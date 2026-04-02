clc; clear; close all;

% Payload (lb)
W_ISR     = 10000;
W_mun     = 20000;
W_crew    = 800;
W_payload = W_ISR + W_mun;

% Mission
Range_leg_nmi = 2000;     % nmi one-way
E_loiter_hr   = 4.0;      % hr loiter on-station
E_reserve_hr  = 0.5;      % hr reserve loiter
fuel_cont     = 0.05;     % contingency / trapped / unusable fraction

% Cruise condition
Mach_cruise = 0.75;
h_cruise_ft = 35000;

% Aero assumptions
LD_cruise = 20;
LD_loiter = 20;

% TSFC in 1/hr for Breguet
c_cruise = 0.55;
c_loiter = 0.50;

% Field length requirements (ft)
TOFL_req_ft = 6000;
LFL_req_ft  = 5000;

% High-lift assumptions
CLmax_TO = 2.2;
CLmax_L  = 2.0;

% Landing mass fraction
mML_over_mMTO = 0.78;

% Runway density ratios sigma = rho/rhoSL
sigma_TO = 1.0;
sigma_L  = 1.0;

% Raymer empty weight regression
A = 0.93;
C = -0.07;

% Engine selection
LEAP_rating_lbf_SL = 30000;   % per engine SL static
nEng               = 2;

% Fixed selected wing area
S_fixed_ft2 = 2250;
S_fixed_m2  = S_fixed_ft2 * 0.09290304;

% Sizing margins
thrust_margin = 1.20;

% Supercritical sizing intent
CL_design_SC2 = 0.40;

% Outer-loop initial guess + tolerances
W0_guess = 170000;
tol      = 1e-4;

% Turn-ons
do_cruise_thrust_check = true;
do_constraint_plot     = true;
do_ceiling_check       = true;

% Ceiling settings
use_cruise_weight_for_ceiling = false;   % more conservative
M_ceiling = Mach_cruise;
h_ceiling_max_ft = 60000;
h_ceiling_step_ft = 500;

% Drag polar for ceiling analysis
CD0_clean = 0.022;
e_clean   = 0.78;
AR_wing   = 11.5;

rhoSL     = 1.225;
kt_per_ms = 1.94384;
g0        = 9.80665;

% Loftin/Raymer-style constants
kL  = 0.107;
kTO = 2.34;


for outer = 1:40

    % 1) Mission fuel fraction
    V_cruise_kt = mach_to_tas_knots(Mach_cruise, h_cruise_ft);
    t_leg_hr    = Range_leg_nmi / V_cruise_kt;

    W1_W0 = 0.970;
    W2_W1 = 0.985;
    W3_W2 = 0.995;
    W7_W6 = 0.995;
    W8_W7 = 0.980;

    W4_W3 = exp(-(t_leg_hr * c_cruise) / LD_cruise);
    W5_W4 = exp(-(E_loiter_hr * c_loiter) / LD_loiter);
    W6_W5 = exp(-(t_leg_hr * c_cruise) / LD_cruise);
    W9_W8 = exp(-(E_reserve_hr * c_loiter) / LD_loiter);

    Mff_nom   = W1_W0 * W2_W1 * W3_W2 * W4_W3 * W5_W4 * W6_W5 * W7_W6 * W8_W7 * W9_W8;
    Wf_W0_nom = 1 - Mff_nom;
    Wf_W0     = 1 - (1 - Wf_W0_nom) * (1 - fuel_cont);

    % Mid-cruise weight fraction
    W3_W0 = W1_W0 * W2_W1 * W3_W2;
    WmidCruise_over_W0 = W3_W0 * sqrt(W4_W3);

    % 2) Cruise atmosphere + dynamic pressure
    sigma_cruise = density_ratio_ISA(h_cruise_ft);
    rho_cruise   = rhoSL * sigma_cruise;

    V_cruise_ms  = V_cruise_kt / kt_per_ms;
    q_cruise     = 0.5 * rho_cruise * V_cruise_ms^2;

    % 3) Solve W0 first
    W0   = solve_W0_from_fractions(W_payload + W_crew, Wf_W0, A, C, W0_guess, tol);
    W0_N = lb_to_N(W0);

    % 4) Fixed wing area -> actual wing loading
    % per Raymer Table 5.5 (pg 84) W/S around 120 expected for bomber
    % chapter 19 discuss optimize w/s and t/w
    S_m2  = S_fixed_m2;
    S_ft2 = S_fixed_ft2;

    WS_MTO_Nm2   = W0_N / S_m2;
    WS_MTO_lbft2 = Nm2_to_lbft2(WS_MTO_Nm2);
    mS_MTO       = WS_MTO_Nm2 / g0;

    % 5) Constraint caps for checking only
    WS_MTO_max_L_Nm2 = WS_from_landing_Loftin_SI(LFL_req_ft, sigma_L, CLmax_L, mML_over_mMTO, kL);
    WS_MTO_max_C_Nm2 = (q_cruise * CL_design_SC2) / max(WmidCruise_over_W0, 1e-6);

    WS_MTO_max_L_lbft2 = Nm2_to_lbft2(WS_MTO_max_L_Nm2);
    WS_MTO_max_C_lbft2 = Nm2_to_lbft2(WS_MTO_max_C_Nm2);

    landing_ok  = WS_MTO_Nm2 <= WS_MTO_max_L_Nm2;
    cruiseCL_ok = WS_MTO_Nm2 <= WS_MTO_max_C_Nm2;

    % 6) Takeoff field length -> required T/W
    TW_req    = TW_from_takeoff_Loftin_SI(TOFL_req_ft, sigma_TO, CLmax_TO, mS_MTO, kTO);
    TW_design = thrust_margin * TW_req;

    % 7) Thrust
    T_total_lbf = TW_design * W0;
    T_each_lbf  = T_total_lbf / nEng;

    outer_err = abs(W0 - W0_guess) / max(W0_guess,1);
    W0_guess  = W0;

    if outer_err < 5e-4
        break;
    end
end

We = (A * W0^C) * W0;
Wf = Wf_W0 * W0;

engine_ok_SL = (T_each_lbf <= LEAP_rating_lbf_SL);

%% CRUISE THRUST CHECK

if do_cruise_thrust_check
    W_cruise_lbf     = WmidCruise_over_W0 * W0;
    D_req_total_lbf  = W_cruise_lbf / LD_cruise;
    T_req_each_lbf   = D_req_total_lbf / nEng;

    lapse_cruise     = thrust_lapse_simple(h_cruise_ft, Mach_cruise);
    T_avail_each_lbf = LEAP_rating_lbf_SL * lapse_cruise;

    engine_ok_cruise = (T_avail_each_lbf >= 1.05 * T_req_each_lbf);
else
    lapse_cruise = NaN;
    T_req_each_lbf = NaN;
    T_avail_each_lbf = NaN;
    engine_ok_cruise = true;
end

%% SPEEDS

rho_TO = rhoSL * sigma_TO;
rho_L  = rhoSL * sigma_L;

Vs_L_ms  = sqrt(2 * WS_MTO_Nm2 * mML_over_mMTO / (rho_L * CLmax_L));
Vapp_ms  = 1.30 * Vs_L_ms;

Vs_TO_ms = sqrt(2 * WS_MTO_Nm2 / (rho_TO * CLmax_TO));
V2_ms    = 1.20 * Vs_TO_ms;

Vs_L_kt  = Vs_L_ms  * kt_per_ms;
Vapp_kt  = Vapp_ms  * kt_per_ms;
Vs_TO_kt = Vs_TO_ms * kt_per_ms;
V2_kt    = V2_ms    * kt_per_ms;

%% CRUISE CONSISTENCY

W_cruise_N     = lb_to_N(WmidCruise_over_W0 * W0);
WS_cruise_Nm2  = W_cruise_N / S_m2;
CL_cruise_req  = WS_cruise_Nm2 / q_cruise;
CD_cruise_req  = CL_cruise_req / LD_cruise;

%% RESULTS

fprintf('\n---- Mission ----\n');
fprintf('Range leg (one-way)         : %.0f nmi\n', Range_leg_nmi);
fprintf('Loiter / Reserve            : %.2f / %.2f hr\n', E_loiter_hr, E_reserve_hr);
fprintf('TSFC cruise / loiter        : %.2f / %.2f 1/hr\n', c_cruise, c_loiter);
fprintf('L/D cruise / loiter         : %.1f / %.1f\n', LD_cruise, LD_loiter);
fprintf('Fuel contingency            : %.1f %%\n', 100 * fuel_cont);

fprintf('\n---- Weights ----\n');
fprintf('Fuel fraction Wf/W0         : %.3f\n', Wf_W0);
fprintf('W0 (MTOW)                   : %.0f lb\n', W0);
fprintf('We (empty)                  : %.0f lb\n', We);
fprintf('Wf (fuel)                   : %.0f lb\n', Wf);
fprintf('Payload + crew              : %.0f lb\n', W_payload + W_crew);

fprintf('\n---- Fixed Wing Area ----\n');
fprintf('Selected wing area          : %.0f ft^2\n', S_ft2);
fprintf('Actual W/S at MTO           : %.1f lb/ft^2\n', WS_MTO_lbft2);
fprintf('Landing W/S cap             : %.1f lb/ft^2\n', WS_MTO_max_L_lbft2);
fprintf('Cruise-CL W/S cap           : %.1f lb/ft^2\n', WS_MTO_max_C_lbft2);
fprintf('Landing check               : %s\n', ternary(landing_ok,'OK','FAIL'));
fprintf('Cruise CL check             : %s\n', ternary(cruiseCL_ok,'OK','FAIL'));

fprintf('\n---- Thrust Sizing ----\n');
fprintf('(T/W)_req from TOFL         : %.3f\n', TW_req);
fprintf('Thrust margin               : %.2f\n', thrust_margin);
fprintf('(T/W)_design                : %.3f\n', TW_design);
fprintf('Total thrust required       : %.0f lbf\n', T_total_lbf);
fprintf('Thrust per engine (n=%d)    : %.0f lbf\n', nEng, T_each_lbf);

fprintf('\n---- Engine checks ----\n');
fprintf('LEAP rating (SL static)     : %.0f lbf per engine\n', LEAP_rating_lbf_SL);
fprintf('SL rating check             : %s\n', ternary(engine_ok_SL,'OK','FAIL'));

if do_cruise_thrust_check
    fprintf('Cruise lapse (simple)       : %.3f\n', lapse_cruise);
    fprintf('Req thrust/engine (cruise)  : %.0f lbf\n', T_req_each_lbf);
    fprintf('Avail thrust/engine (cruise): %.0f lbf\n', T_avail_each_lbf);
    fprintf('Cruise thrust check         : %s\n', ternary(engine_ok_cruise,'OK','FAIL'));
end

fprintf('\n---- Speeds (sigma-scaled) ----\n');
fprintf('Landing stall Vs_L          : %.1f kt\n', Vs_L_kt);
fprintf('Approach Vapp=1.3Vs         : %.1f kt\n', Vapp_kt);
fprintf('TO stall Vs_TO              : %.1f kt\n', Vs_TO_kt);
fprintf('TO safety V2=1.2Vs          : %.1f kt\n', V2_kt);

fprintf('\n---- Cruise consistency ----\n');
fprintf('Mid-cruise W/W0             : %.3f\n', WmidCruise_over_W0);
fprintf('Required CL at cruise       : %.3f (target clamp %.2f)\n', CL_cruise_req, CL_design_SC2);
fprintf('Required CD (from L/D=%.1f) : %.4f\n', LD_cruise, CD_cruise_req);
fprintf('Outer iterations            : %d\n', outer);

%% CONSTRAINT DIAGRAM

WS_sweep_lbft2 = linspace(40, 220, 400);
ROC_req_fpm    = 1500;
LD_climb       = 20;
eta_climb      = 0.90;
Vclimb_kt      = max(250, V2_kt);

WS_sweep_Nm2 = WS_sweep_lbft2 * 47.88025898;
mS_sweep     = WS_sweep_Nm2 / g0;

sTO_m = ft_to_m(TOFL_req_ft);
TW_TO = (mS_sweep .* (kTO * sigma_TO)) ./ (sTO_m * CLmax_TO);

TW_cruise_SL = (WmidCruise_over_W0 / LD_cruise) / max(lapse_cruise, 1e-6);

ROC_ms    = ROC_req_fpm * 0.00508;
Vclimb_ms = Vclimb_kt / kt_per_ms;

TW_climb_alt = eta_climb * ((1 / LD_climb) + (ROC_ms / max(Vclimb_ms,1e-6)));
TW_climb_SL  = TW_climb_alt / max(lapse_cruise, 1e-6);

TW_env = max([TW_TO(:), TW_cruise_SL * ones(numel(TW_TO),1), TW_climb_SL * ones(numel(TW_TO),1)], [], 2);

if do_constraint_plot
    figure; hold on; grid on;
    plot(WS_sweep_lbft2, TW_TO, 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_cruise_SL * ones(size(WS_sweep_lbft2)), 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_climb_SL * ones(size(WS_sweep_lbft2)), 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_env, 'k--', 'LineWidth', 2.0);

    xline(WS_MTO_max_L_lbft2, '-.', 'LineWidth', 1.5);
    xline(WS_MTO_max_C_lbft2, ':', 'LineWidth', 1.5);
    xline(WS_MTO_lbft2, '--', 'LineWidth', 1.5);

    plot(WS_MTO_lbft2, TW_design, 'o', 'MarkerSize', 8, 'LineWidth', 1.8);

    xlabel('W/S at MTO (lb/ft^2)');
    ylabel('Required T/W at SL (static equivalent)');
    title('Constraint Diagram with Fixed Wing Area');
    legend({'Takeoff (TOFL)', 'Cruise', sprintf('Climb (ROC=%d fpm)',ROC_req_fpm), ...
            'Envelope (max)', 'Landing W/S cap', 'Cruise CL cap', 'Actual fixed-area W/S', 'Design point'}, ...
            'Location','best');
end

%% CEILING ANALYSIS 

if do_ceiling_check

    if use_cruise_weight_for_ceiling
        W_ceiling_lbf = WmidCruise_over_W0 * W0;
        ceiling_weight_label = 'Cruise weight';
    else
        W_ceiling_lbf = W0;
        ceiling_weight_label = 'MTOW';
    end

    h_grid_ft = 0:h_ceiling_step_ft:h_ceiling_max_ft;

    T_avail_total_lbf = zeros(size(h_grid_ft));
    T_req_total_lbf   = zeros(size(h_grid_ft));
    ROC_fpm           = zeros(size(h_grid_ft));
    CL_grid           = zeros(size(h_grid_ft));
    CD_grid           = zeros(size(h_grid_ft));
    LD_grid           = zeros(size(h_grid_ft));

    k_induced = 1 / (pi * e_clean * AR_wing);
    W_ceiling_N = lb_to_N(W_ceiling_lbf);

    for i = 1:length(h_grid_ft)
        h_ft = h_grid_ft(i);
        sigma_i = density_ratio_ISA(h_ft);
        rho_i   = rhoSL * sigma_i;

        V_kt = mach_to_tas_knots(M_ceiling, h_ft);
        V_ms = V_kt / kt_per_ms;

        q_i = 0.5 * rho_i * V_ms^2;

        CL_grid(i) = W_ceiling_N / (q_i * S_m2);
        CD_grid(i) = CD0_clean + k_induced * CL_grid(i)^2;
        LD_grid(i) = CL_grid(i) / CD_grid(i);

        T_req_total_lbf(i)   = W_ceiling_lbf / LD_grid(i);
        T_avail_total_lbf(i) = nEng * LEAP_rating_lbf_SL * thrust_lapse_simple(h_ft, M_ceiling);

        ROC_fpm(i) = ((T_avail_total_lbf(i) - T_req_total_lbf(i)) / W_ceiling_lbf) * V_ms * 196.8504;
    end

    diffT = T_avail_total_lbf - T_req_total_lbf;
    idx_abs = find(diffT(1:end-1) >= 0 & diffT(2:end) <= 0, 1, 'first');

    if ~isempty(idx_abs)
        h1 = h_grid_ft(idx_abs);
        h2 = h_grid_ft(idx_abs+1);
        d1 = diffT(idx_abs);
        d2 = diffT(idx_abs+1);
        absolute_ceiling_ft = h1 + (0 - d1) * (h2 - h1) / (d2 - d1);
    else
        absolute_ceiling_ft = NaN;
    end

    ROC_service_criterion = 100;
    idx_serv = find(ROC_fpm(1:end-1) >= ROC_service_criterion & ROC_fpm(2:end) <= ROC_service_criterion, 1, 'first');

    if ~isempty(idx_serv)
        h1 = h_grid_ft(idx_serv);
        h2 = h_grid_ft(idx_serv+1);
        r1 = ROC_fpm(idx_serv);
        r2 = ROC_fpm(idx_serv+1);
        service_ceiling_ft = h1 + (ROC_service_criterion - r1) * (h2 - h1) / (r2 - r1);
    else
        service_ceiling_ft = NaN;
    end

    fprintf('\n---- Ceiling Analysis ----\n');
    fprintf('Ceiling weight basis        : %s\n', ceiling_weight_label);
    fprintf('Ceiling analysis weight     : %.0f lb\n', W_ceiling_lbf);
    fprintf('Ceiling analysis Mach       : %.2f\n', M_ceiling);
    fprintf('CD0 used                    : %.4f\n', CD0_clean);
    fprintf('Oswald e used               : %.2f\n', e_clean);
    fprintf('Aspect ratio used           : %.2f\n', AR_wing);

    if ~isnan(service_ceiling_ft)
        fprintf('Service Ceiling             : %.2f ft\n', service_ceiling_ft);
    else
        fprintf('Service Ceiling             : Not reached in sweep\n');
    end

    if ~isnan(absolute_ceiling_ft)
        fprintf('Absolute Ceiling            : %.2f ft\n', absolute_ceiling_ft);
    else
        fprintf('Absolute Ceiling            : Not reached in sweep\n');
    end

    fprintf('30,000 ft requirement       : %s\n', ternary(~isnan(service_ceiling_ft) && service_ceiling_ft >= 30000,'PASS','FAIL'));

    figure; hold on; grid on;
    plot(h_grid_ft, T_avail_total_lbf, 'LineWidth', 1.8);
    plot(h_grid_ft, T_req_total_lbf, 'LineWidth', 1.8);
    xlabel('Altitude (ft)');
    ylabel('Thrust (lbf)');
    title('Available and Required Thrust vs Altitude');
    legend('Available thrust','Required thrust','Location','best');

    figure; hold on; grid on;
    plot(h_grid_ft, ROC_fpm, 'LineWidth', 1.8);
    yline(ROC_service_criterion, '--', 'LineWidth', 1.5);
    yline(300, ':', 'LineWidth', 1.5);
    xlabel('Altitude (ft)');
    ylabel('Rate of Climb (ft/min)');
    title('Rate of Climb vs Altitude');
    legend('ROC',sprintf('%d ft/min criterion',ROC_service_criterion),'300 ft/min','Location','best');
end

%% FUNCTIONS

function W0 = solve_W0_from_fractions(W_fixed, Wf_W0, A, C, W0_init, tol)
    W0  = W0_init;
    err = 1;
    it  = 0;
    while err > tol && it < 800
        We_W0 = A * W0^C;
        denom = 1 - We_W0 - Wf_W0;
        if denom <= 0.02
            error('Infeasible fractions: 1 - We/W0 - Wf/W0 = %.4f.', denom);
        end
        W0_new = W_fixed / denom;
        err    = abs(W0_new - W0) / max(W0,1);
        W0     = W0_new;
        it     = it + 1;
    end
end

function [WS_lim_Nm2, Sref_m2, details] = wingSizingRaymerConstraintSI( ...
    TOFL_ft, LFL_ft, ...
    sigma_TO, sigma_L, ...
    CLmax_TO, CLg_TO, CD0_TO, K_TO, mu_TO, TW_TO, ...
    CLmax_L, mL_over_mTO, Sa_ft, landingFactor, WTO_N)
%WINGSIZINGRAYMERCONSTRAINTSI
% Sizes wing loading and wing area from Raymer-style takeoff and landing constraints.
%
% IMPORTANT:
%   - Takeoff uses the Raymer ground-roll equation only.
%   - If TOFL_ft is a full field length to clear a 50 ft obstacle, you should
%     add rotation + transition + climb separately.
%   - Landing uses Raymer's simplified wing-loading relation.

    if nargin < 14 || isempty(landingFactor)
        landingFactor = 1.0;   % 1.0 = basic Raymer relation
    end

    rhoSL = 1.225; % kg/m^3
    g0    = 9.80665;

    % Takeoff wing-loading limit from Raymer ground-roll equation
    WS_takeoff_max_Nm2 = WS_from_takeoff_Raymer_SI( ...
        TOFL_ft, sigma_TO, CLmax_TO, CLg_TO, CD0_TO, K_TO, mu_TO, TW_TO);

    % Landing wing-loading limit from Raymer approximate landing relation
    WS_landing_max_Nm2 = WS_from_landing_Raymer_SI( ...
        LFL_ft, sigma_L, CLmax_L, mL_over_mTO, Sa_ft, landingFactor);

    % Limiting wing loading
    WS_lim_Nm2 = min(WS_takeoff_max_Nm2, WS_landing_max_Nm2);

    % Wing area from takeoff weight
    Sref_m2 = WTO_N / WS_lim_Nm2;

    % Pack outputs
    details = struct();
    details.WS_takeoff_max_Nm2 = WS_takeoff_max_Nm2;
    details.WS_landing_max_Nm2 = WS_landing_max_Nm2;
    details.limiting_case = "takeoff";
    if WS_landing_max_Nm2 < WS_takeoff_max_Nm2
        details.limiting_case = "landing";
    end
    details.WTO_N = WTO_N;
    details.Sref_m2 = Sref_m2;
end

function WS_max_Nm2 = WS_from_takeoff_Raymer_SI( ...
    TOFL_ft, sigma_TO, CLmax_TO, CLg_TO, CD0_TO, K_TO, mu_TO, TW_TO)
% Solve for max W/S (N/m^2) such that Raymer takeoff ground roll <= TOFL_ft.
%
% Uses:
%   Sg = 1/(2 g Ka) * ln(1 + (Ka/Kt) * Vto^2)
%   Kt = T/W - mu
%   Ka = rho/(2*(W/S)) * (mu*CLg - CD0 - K*CLg^2)
%   Vto = 1.1 * sqrt(2*(W/S)/(rho*CLmax_TO))

    g0    = 9.80665;
    rhoSL = 1.225;           % kg/m^3
    rho   = sigma_TO * rhoSL; % takeoff density

    target_m = ft_to_m(TOFL_ft);

    if TW_TO <= mu_TO
        error('Takeoff infeasible: T/W must exceed mu_TO.');
    end

    % Residual function: Sg(WS) - target
    f = @(WS) takeoff_residual(WS, target_m, rho, CLmax_TO, CLg_TO, ...
                              CD0_TO, K_TO, mu_TO, TW_TO, g0);

    % Bracket the solution
    WS_lo = 1.0;      % N/m^2
    WS_hi = 2.0e5;    % N/m^2

    f_lo = f(WS_lo);
    f_hi = f(WS_hi);

    % Expand upper bound until we bracket the root
    nGrow = 0;
    while ~(isfinite(f_lo) && isfinite(f_hi) && f_lo <= 0 && f_hi >= 0)
        if ~isfinite(f_hi) || f_hi < 0
            WS_hi = WS_hi * 2.0;
            f_hi = f(WS_hi);
        elseif f_lo > 0
            WS_lo = WS_lo / 2.0;
            f_lo = f(WS_lo);
        end

        nGrow = nGrow + 1;
        if nGrow > 80 || WS_hi > 5.0e7 || WS_lo < 1.0e-8
            error('Could not bracket the takeoff wing-loading solution.');
        end
    end

    % Bisection
    for k = 1:100
        WS_mid = 0.5 * (WS_lo + WS_hi);
        f_mid  = f(WS_mid);

        if ~isfinite(f_mid)
            WS_lo = WS_mid;
            continue;
        end

        if f_mid > 0
            WS_hi = WS_mid;
        else
            WS_lo = WS_mid;
        end
    end

    WS_max_Nm2 = 0.5 * (WS_lo + WS_hi);
end

function r = takeoff_residual(WS, target_m, rho, CLmax_TO, CLg_TO, ...
                              CD0_TO, K_TO, mu_TO, TW_TO, g0)
% Residual for root finding: Sg(WS) - target

    if WS <= 0
        r = Inf;
        return;
    end

    Kt = TW_TO - mu_TO;
    if Kt <= 0
        r = Inf;
        return;
    end

    Vto = 1.1 * sqrt(2.0 * WS / (rho * CLmax_TO));

    Ka = (rho / (2.0 * WS)) * (mu_TO * CLg_TO - CD0_TO - K_TO * CLg_TO^2);

    % Fallback if Ka is numerically tiny
    if abs(Ka) < 1e-12
        Sg = Vto^2 / (2.0 * g0 * Kt);
        r = Sg - target_m;
        return;
    end

    arg = 1.0 + (Ka / Kt) * Vto^2;

    if arg <= 0
        r = Inf;
        return;
    end

    Sg = log(arg) / (2.0 * g0 * Ka);
    r  = Sg - target_m;
end

function WS_TO_max_Nm2 = WS_from_landing_Raymer_SI( ...
    LFL_ft, sigma_L, CLmax_L, mL_over_mTO, Sa_ft, landingFactor)
% Raymer landing approximation:
%   S_landing = landingFactor * 80 * (W/S)/(sigma * CLmax) + Sa
%
% This is the classic constraint-analysis form in English units.
% It returns W/S in N/m^2 after converting from lb/ft^2.

    if mL_over_mTO <= 0 || mL_over_mTO > 1.0
        error('mL_over_mTO must be in (0, 1].');
    end

    if LFL_ft <= Sa_ft
        error('Landing field length must be greater than Sa_ft.');
    end

    % Solve in lb/ft^2 first because the Raymer constant "80" is in English units.
    WS_L_lbft2 = (LFL_ft - Sa_ft) * sigma_L * CLmax_L / (80.0 * landingFactor);

    % Convert to N/m^2
    WS_L_Nm2 = WS_L_lbft2 * 47.88025898;

    % Convert from landing weight to takeoff weight
    WS_TO_max_Nm2 = WS_L_Nm2 / mL_over_mTO;
end

function m = ft_to_m(ft)
    m = ft * 0.3048;
end

function V = mach_to_tas_knots(M, h_ft)
    gamma = 1.4;
    R = 287.05;
    T0 = 288.15;
    aL = -0.0065;
    h  = h_ft * 0.3048;
    if h <= 11000
        T = T0 + aL * h;
    else
        T11 = T0 + aL * 11000;
        T = T11;
    end
    a    = sqrt(gamma * R * T);
    V_ms = M * a;
    V    = V_ms * 1.94384;
end

function lapse = thrust_lapse_simple(h_ft, M)
    sigma = density_ratio_ISA(h_ft);
    lapse = (sigma^1.0) * (1 - 0.35 * M);
    lapse = max(min(lapse, 1.0), 0.05);
end

function sigma = density_ratio_ISA(h_ft)
    rho0 = 1.225;
    T0   = 288.15;
    p0   = 101325;
    g0   = 9.80665;
    R    = 287.05;
    aL   = -0.0065;
    h = h_ft * 0.3048;
    if h <= 11000
        T = T0 + aL * h;
        p = p0 * (T / T0)^(-g0 / (aL * R));
    else
        T11 = T0 + aL * 11000;
        p11 = p0 * (T11 / T0)^(-g0 / (aL * R));
        T = T11;
        p = p11 * exp(-g0 * (h - 11000) / (R * T));
    end
    rho   = p / (R * T);
    sigma = rho / rho0;
end

function m = ft_to_m(ft)
    m = ft * 0.3048;
end

function N = lb_to_N(lb)
    N = lb * 4.4482216152605;
end

function lbft2 = Nm2_to_lbft2(Nm2)
    lbft2 = Nm2 / 47.88025898;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end