
clc; clear; close all;

% Payload (lb)
W_ISR     = 10000;
W_mun     = 20000;
W_crew    = 800;
W_payload = W_ISR + W_mun;

% Mission
Range_leg_nmi = 2000;     % nmi one-way
E_loiter_hr   = 4.0;      % hr loiter on-station
E_reserve_hr  = 0.75;     % hr reserve loiter
fuel_cont     = 0.05;     % contingency/trapped/unusable fraction

% Cruise condition
Mach_cruise = 0.75;
h_cruise_ft = 30000;

% Aero assumptions
LD_cruise = 18;
LD_loiter = 16;

% TSFC in 1/hr for Breguet
c_cruise = 0.55;
c_loiter = 0.50;

% Field length requirements (ft)
TOFL_req_ft = 6000;
LFL_req_ft  = 5000;

% High-lift assumptions
CLmax_TO = 2.5;
CLmax_L  = 2.8;

% Landing mass fraction
mML_over_mMTO = 0.78;

% Runway density ratios sigma = rho/rhoSL (set <1 for hot/high)
sigma_TO = 1.0;
sigma_L  = 1.0;

% Raymer empty weight regression: We/W0 = A * W0^C  (W0 in lb)
A = 0.93;
C = -0.07;

% Engine selection (conceptual)
LEAP_rating_lbf_SL = 31000;  % per engine SL static
nEng               = 2;

% Sizing margins
thrust_margin = 1.20;  % multiplies TOFL-based T/W
margin_WS     = 0.95;  % keeps away from W/S caps

% Supercritical sizing intent (cruise CL clamp)
CL_design_SC2 = 0.40;

% Outer-loop initial guess + tolerances
W0_guess = 170000;     % lb
tol      = 1e-4;

% Turn-ons
do_cruise_thrust_check = true;
do_constraint_plot     = true;

rhoSL     = 1.225;       % kg/m^3
kt_per_ms = 1.94384;
g0        = 9.80665;     % m/s^2

% Loftin/Raymer-style constants used in your prior code (kept isolated!)
kL  = 0.107;     % landing constant for SI form used here
kTO = 2.34;      % takeoff constant for SI form used here

for outer = 1:40

    % 1) Mission fuel fraction (segment + Breguet)
    V_cruise_kt = mach_to_tas_knots(Mach_cruise, h_cruise_ft);
    t_leg_hr    = Range_leg_nmi / V_cruise_kt;

    % Segment fractions (adjust if course uses different)
    W1_W0 = 0.970;   % warmup/taxi
    W2_W1 = 0.985;   % takeoff
    W3_W2 = 0.995;   % climb
    W7_W6 = 0.995;   % descent
    W8_W7 = 0.980;   % landing

    % Breguet segments
    W4_W3 = exp(-(t_leg_hr * c_cruise) / LD_cruise);      % cruise out
    W5_W4 = exp(-(E_loiter_hr * c_loiter) / LD_loiter);   % loiter
    W6_W5 = exp(-(t_leg_hr * c_cruise) / LD_cruise);      % cruise back
    W9_W8 = exp(-(E_reserve_hr * c_loiter) / LD_loiter);  % reserve loiter

    Mff_nom   = W1_W0 * W2_W1 * W3_W2 * W4_W3 * W5_W4 * W6_W5 * W7_W6 * W8_W7 * W9_W8;
    Wf_W0_nom = 1 - Mff_nom;
    Wf_W0     = 1 - (1 - Wf_W0_nom) * (1 - fuel_cont);   % add trapped/contingency

    % Consistent mid-cruise weight fraction (mid outbound cruise)
    W3_W0 = W1_W0 * W2_W1 * W3_W2;
    WmidCruise_over_W0 = W3_W0 * sqrt(W4_W3);

    %2) Cruise atmosphere + dynamic pressure
    sigma_cruise = density_ratio_ISA(h_cruise_ft);
    rho_cruise   = rhoSL * sigma_cruise;

    V_cruise_ms  = V_cruise_kt / kt_per_ms;
    q_cruise     = 0.5 * rho_cruise * V_cruise_ms^2;  % N/m^2

    % 3) Landing field-length constraint -> max W/S at MTO
    WS_MTO_max_L_Nm2 = WS_from_landing_Loftin_SI(LFL_req_ft, sigma_L, CLmax_L, mML_over_mMTO, kL);

    % 4) Cruise CL clamp -> max W/S at MTO
    WS_MTO_max_C_Nm2 = (q_cruise * CL_design_SC2) / max(WmidCruise_over_W0, 1e-6);

    %5) Choose W/S as most restrictive (with margin)
    [WS_MTO_max_Nm2, WS_driver] = min_with_label(WS_MTO_max_L_Nm2, "landing", WS_MTO_max_C_Nm2, "cruise_CL");
    WS_MTO_Nm2   = margin_WS * WS_MTO_max_Nm2;
    WS_MTO_lbft2 = Nm2_to_lbft2(WS_MTO_Nm2);
    mS_MTO       = WS_MTO_Nm2 / g0;   % kg/m^2

    %6) Takeoff field-length constraint -> required T/W at SL (static equivalent)
    TW_req    = TW_from_takeoff_Loftin_SI(TOFL_req_ft, sigma_TO, CLmax_TO, mS_MTO, kTO);
    TW_design = thrust_margin * TW_req;

    %7) Solve W0 from payload+crew, fuel fraction, empty fraction regression
    W0 = solve_W0_from_fractions(W_payload + W_crew, Wf_W0, A, C, W0_guess, tol);

    %8) Geometry + thrust from W/S and T/W
    W0_N  = lb_to_N(W0);
    S_m2  = W0_N / WS_MTO_Nm2;
    S_ft2 = m2_to_ft2(S_m2);

    T_total_lbf = TW_design * W0;
    T_each_lbf  = T_total_lbf / nEng;

    %9) Outer convergence
    outer_err = abs(W0 - W0_guess) / max(W0_guess,1);
    W0_guess  = W0;

    if outer_err < 5e-4
        break;
    end
end

We = (A * W0^C) * W0;
Wf = Wf_W0 * W0;

engine_ok_SL = (T_each_lbf <= LEAP_rating_lbf_SL);

% Cruise thrust sanity check (rough): D ≈ W/(L/D)
if do_cruise_thrust_check
    W_cruise_lbf     = WmidCruise_over_W0 * W0;
    D_req_total_lbf  = W_cruise_lbf / LD_cruise;
    T_req_each_lbf   = D_req_total_lbf / nEng;

    lapse_cruise     = thrust_lapse_simple(h_cruise_ft, Mach_cruise);
    T_avail_each_lbf = LEAP_rating_lbf_SL * lapse_cruise;

    engine_ok_cruise = (T_avail_each_lbf >= 1.05*T_req_each_lbf);
else
    lapse_cruise = NaN; T_req_each_lbf = NaN; T_avail_each_lbf = NaN; engine_ok_cruise = true;
end

% Runway densities
rho_TO = rhoSL * sigma_TO;
rho_L  = rhoSL * sigma_L;

% Landing speeds at ML (use MTO W/S scaled to ML by mML/mMTO)
Vs_L_ms  = sqrt( 2*WS_MTO_Nm2 * mML_over_mMTO / (rho_L * CLmax_L) );
Vapp_ms  = 1.30 * Vs_L_ms;

% Takeoff speeds at MTO
Vs_TO_ms = sqrt( 2*WS_MTO_Nm2 / (rho_TO * CLmax_TO) );
V2_ms    = 1.20 * Vs_TO_ms;

Vs_L_kt  = Vs_L_ms  * kt_per_ms;
Vapp_kt  = Vapp_ms  * kt_per_ms;
Vs_TO_kt = Vs_TO_ms * kt_per_ms;
V2_kt    = V2_ms    * kt_per_ms;

% Cruise CL and CD required from L/D
W_cruise_N     = lb_to_N(WmidCruise_over_W0 * W0);
WS_cruise_Nm2  = W_cruise_N / S_m2;
CL_cruise_req  = WS_cruise_Nm2 / q_cruise;
CD_cruise_req  = CL_cruise_req / LD_cruise;

fprintf('\n---- Mission ----\n');
fprintf('Range leg (one-way)         : %.0f nmi\n', Range_leg_nmi);
fprintf('Loiter / Reserve            : %.2f / %.2f hr\n', E_loiter_hr, E_reserve_hr);
fprintf('TSFC cruise / loiter        : %.2f / %.2f 1/hr\n', c_cruise, c_loiter);
fprintf('L/D cruise / loiter         : %.1f / %.1f\n', LD_cruise, LD_loiter);
fprintf('Fuel contingency            : %.1f %%\n', 100*fuel_cont);

fprintf('\n---- Weights ----\n');
fprintf('Fuel fraction Wf/W0         : %.3f\n', Wf_W0);
fprintf('W0 (MTOW)                   : %.0f lb\n', W0);
fprintf('We (empty)                  : %.0f lb\n', We);
fprintf('Wf (fuel)                   : %.0f lb\n', Wf);
fprintf('Payload + crew              : %.0f lb\n', W_payload + W_crew);

fprintf('\n---- Wing Loading / Area ----\n');
fprintf('W/S (MTO)                   : %.1f lb/ft^2  (driver: %s)\n', WS_MTO_lbft2, WS_driver);
fprintf('Wing area S                 : %.0f ft^2\n', S_ft2);

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

if ~engine_ok_SL || (do_cruise_thrust_check && ~engine_ok_cruise)
    fprintf('\n*** DESIGN NOTES ***\n');
    fprintf('- If thrust too high: lower W/S (bigger wing), increase CLmax_TO, relax TOFL, increase engine count, or raise L/D.\n');
    fprintf('- If cruise thrust fails: increase lapse fidelity, reduce drag assumptions, or increase thrust margin.\n');
end


% W/S caps in lb/ft^2 for reporting/plot
WS_MTO_max_L_lbft2 = Nm2_to_lbft2(WS_MTO_max_L_Nm2);
WS_MTO_max_C_lbft2 = Nm2_to_lbft2(WS_MTO_max_C_Nm2);

% ---- Mission segment table ----
seg_names = {'Warmup/Taxi','Takeoff','Climb','Cruise Out','Loiter','Cruise Back','Descent','Landing','Reserve Loiter'};
seg_frac  = [W1_W0, W2_W1, W3_W2, W4_W3, W5_W4, W6_W5, W7_W6, W8_W7, W9_W8];

cum = ones(size(seg_frac));
cum(1) = seg_frac(1);
for i = 2:numel(seg_frac)
    cum(i) = cum(i-1) * seg_frac(i);
end
seg_fuel = 1 - seg_frac;

T_mission = table(seg_names(:), seg_frac(:), cum(:), seg_fuel(:), ...
    'VariableNames', {'Segment','Wi_over_Wprev','W_over_W0_cum','FuelFrac_inSegment'});

disp(' ');
disp('================ MISSION SEGMENTS (FRACTIONS) ================');
disp(T_mission);

% ---- Planform estimates from AR sweep ----
AR_list = [8 9 10 11];
disp(' ');
disp('================ PLANFORM GUESSES (from S & AR) ================');
for AR = AR_list
    b_m  = sqrt(AR * S_m2);
    c_m  = S_m2 / b_m;
    fprintf('AR=%4.1f -> span b=%.1f ft | mean chord c=%.1f ft\n', AR, b_m/0.3048, c_m/0.3048);
end

% ---- Fuel volume sanity check (rough) ----
rho_fuel_lb_per_gal = 6.7;   % Jet-A approx 
t_c_guess           = 0.12;  % thickness ratio guess
tankable_fraction   = 0.55;  % usable tank fraction

Vfuel_gal = Wf / rho_fuel_lb_per_gal;

AR_ref = 9;
b_ref  = sqrt(AR_ref * S_m2);
cbar   = S_m2 / b_ref;
t_avg  = t_c_guess * cbar;

Vwing_m3_solid = S_m2 * t_avg;                 % crude "solid" volume
Vtank_m3       = tankable_fraction * Vwing_m3_solid;
Vtank_gal      = Vtank_m3 * 264.172;

disp(' ');
disp('================ FUEL VOLUME SANITY (rough) ==================');
fprintf('Fuel weight Wf              : %.0f lb\n', Wf);
fprintf('Fuel volume (Jet-A)         : %.0f gal (%.1f lb/gal)\n', Vfuel_gal, rho_fuel_lb_per_gal);
fprintf('Wing tank estimate          : %.0f gal (AR=%.1f, t/c=%.2f, tankable=%.2f)\n', ...
    Vtank_gal, AR_ref, t_c_guess, tankable_fraction);

if Vtank_gal < 0.9*Vfuel_gal
    fprintf('NOTE: Fuel volume looks TIGHT → consider thicker wing, body tanks, or more tankable fraction.\n');
elseif Vtank_gal > 1.3*Vfuel_gal
    fprintf('NOTE: Fuel volume looks OK (or generous) under these assumptions.\n');
else
    fprintf('NOTE: Fuel volume looks plausible but sensitive to t/c and tankable fraction.\n');
end

% ---- Constraint diagram ----
% Settings
WS_sweep_lbft2 = linspace(40, 220, 400);
ROC_req_fpm    = 1500;    % conceptual climb requirement
LD_climb       = 16;      % conceptual climb L/D
eta_climb      = 0.90;    % nonideal factor
Vclimb_kt      = max(250, V2_kt);  % simple choice: at least 250 kt or V2

WS_sweep_Nm2 = WS_sweep_lbft2 * 47.88025898;
mS_sweep     = WS_sweep_Nm2 / g0;

% Takeoff curve
sTO_m = ft_to_m(TOFL_req_ft);
TW_TO = (mS_sweep .* (kTO * sigma_TO)) ./ (sTO_m * CLmax_TO);

% Cruise curve (SL static equivalent): T/W_SL ≈ (Wcr/W0)/(L/D)/lapse
TW_cruise_SL = (WmidCruise_over_W0 / LD_cruise) / max(lapse_cruise, 1e-6);

% Climb curve: T/W >= 1/(L/D_climb) + ROC/V
ROC_ms   = ROC_req_fpm * 0.00508;
Vclimb_ms = Vclimb_kt / kt_per_ms;

TW_climb_alt = eta_climb * ( (1/LD_climb) + (ROC_ms / max(Vclimb_ms,1e-6)) );
TW_climb_SL  = TW_climb_alt / max(lapse_cruise, 1e-6);

% Envelope
TW_env = max([TW_TO(:), TW_cruise_SL*ones(numel(TW_TO),1), TW_climb_SL*ones(numel(TW_TO),1)], [], 2);

fprintf('Landing W/S cap             : %.1f lb/ft^2\n', WS_MTO_max_L_lbft2);
fprintf('Cruise-CL W/S cap           : %.1f lb/ft^2\n', WS_MTO_max_C_lbft2);

if do_constraint_plot
    figure; hold on; grid on;
    plot(WS_sweep_lbft2, TW_TO, 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_cruise_SL*ones(size(WS_sweep_lbft2)), 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_climb_SL*ones(size(WS_sweep_lbft2)), 'LineWidth', 1.6);
    plot(WS_sweep_lbft2, TW_env, 'k--', 'LineWidth', 2.0);

    xline(WS_MTO_max_L_lbft2, '-.', 'LineWidth', 1.5);
    xline(WS_MTO_max_C_lbft2, ':', 'LineWidth', 1.5);

    plot(WS_MTO_lbft2, TW_design, 'o', 'MarkerSize', 8, 'LineWidth', 1.8);

    xlabel('W/S at MTO (lb/ft^2)');
    ylabel('Required T/W at SL (static equivalent)');
    title('Constraint Diagram: TOFL, Cruise, Climb + W/S Caps (Landing, Cruise-CL)');
    legend({'Takeoff (TOFL)', 'Cruise', sprintf('Climb (ROC=%d fpm)',ROC_req_fpm), ...
            'Envelope (max)', 'Landing W/S cap', 'Cruise CL cap', 'Design point'}, ...
            'Location','best');
end

%% functs

function W0 = solve_W0_from_fractions(W_fixed, Wf_W0, A, C, W0_init, tol)
    W0  = W0_init;
    err = 1; it = 0;

    while err > tol && it < 800
        We_W0 = A * W0^C;
        denom = 1 - We_W0 - Wf_W0;

        if denom <= 0.02
            error('Infeasible fractions: 1 - We/W0 - Wf/W0 = %.4f. Check A,C and/or mission fuel fraction.', denom);
        end

        W0_new = W_fixed / denom;
        err    = abs(W0_new - W0) / max(W0,1);
        W0     = W0_new;
        it     = it + 1;
    end
end

function WS_MTO_max_Nm2 = WS_from_landing_Loftin_SI(LFL_ft, sigma_L, CLmax_L, mML_over_mMTO, kL)
    g0 = 9.80665;
    sL_m = ft_to_m(LFL_ft);

    % Your chosen Loftin/Raymer SI form (isolated here)
    mS_ML_max      = kL * sigma_L * CLmax_L * sL_m;   % kg/m^2 at ML
    mS_MTO_max_L   = mS_ML_max / mML_over_mMTO;       % kg/m^2 at MTO
    WS_MTO_max_Nm2 = mS_MTO_max_L * g0;               % N/m^2
end

function TW_req = TW_from_takeoff_Loftin_SI(TOFL_ft, sigma_TO, CLmax_TO, mS_MTO, kTO)
    sTO_m = ft_to_m(TOFL_ft);

    % Your chosen Loftin/Raymer SI form (isolated here)
    TW_req = mS_MTO * (kTO * sigma_TO) / (sTO_m * CLmax_TO);
end

function [val, label] = min_with_label(a, la, b, lb)
    if a <= b
        val = a; label = la;
    else
        val = b; label = lb;
    end
end

function V = mach_to_tas_knots(M, h_ft)
    gamma = 1.4;
    R = 287.05;      % J/(kg*K)
    T0 = 288.15;     % K
    aL = -0.0065;    % K/m
    h  = h_ft * 0.3048;

    if h <= 11000
        T = T0 + aL*h;
    else
        T11 = T0 + aL*11000;
        T = T11;
    end

    a    = sqrt(gamma*R*T);
    V_ms = M*a;
    V    = V_ms * 1.94384;
end

function lapse = thrust_lapse_simple(h_ft, M)
    sigma = density_ratio_ISA(h_ft);
    lapse = (sigma^0.8) * (1 - 0.30*M);
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
        T = T0 + aL*h;
        p = p0 * (T/T0)^(-g0/(aL*R));
    else
        T11 = T0 + aL*11000;
        p11 = p0 * (T11/T0)^(-g0/(aL*R));
        T = T11;
        p = p11 * exp(-g0*(h-11000)/(R*T));
    end

    rho   = p/(R*T);
    sigma = rho / rho0;
end

function m = ft_to_m(ft), m = ft * 0.3048; end
function ft2 = m2_to_ft2(m2), ft2 = m2 / (0.3048^2); end

function N = lb_to_N(lb)
    N = lb * 4.4482216152605;
end

function lbft2 = Nm2_to_lbft2(Nm2)
    lbft2 = Nm2 / 47.88025898;  % 1 lb/ft^2 = 47.88025898 N/m^2
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
