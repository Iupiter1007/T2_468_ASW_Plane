% AEE 468 Aircraft Design
clc; clear; close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%% SECTION 1 INPUTS
W_ISR     = 10000;
W_mun     = 20000;
W_crew    = 800;
W_payload = W_ISR + W_mun + W_crew;

Range_leg_nmi = 2000;
E_loiter_hr   = 4.0;
E_reserve_hr  = 0.5;
fuel_cont     = 0.10;
Mach_cruise = 0.75;   h_cruise_ft = 35000;
Mach_loiter = 0.40;   h_loiter_ft = 15000;
LD_cruise = 18;   LD_loiter = 20;
c_cruise  = 0.50; c_loiter  = 0.47;   % TSFC (1/hr)
TOFL_req_ft = 6000;
LFL_req_ft  = 5000;
CLmax_TO       = 2.2;
CLmax_L        = 2.2;
W_land_over_W0 = 0.72;
rhoSL     = 1.225;
kt_per_ms = 1.94384;
CL_TO   = 0.8 * CLmax_TO;
TOP_jet = 155;

LEAP_rating_lbf_SL = 27120;
nEng               = 2;

CD0_clean = 0.020;
e_clean   = 0.82;
AR_wing   = 11.5;
k_induced = 1 / (pi * e_clean * AR_wing);

A_ew = 0.93;   % Empty weight regression (Raymer Table 3.1)
C_ew = -0.07;

W0_guess      = 200000;
tol           = 1e-6;
thrust_margin = 1.10;

CL_design_SC = 0.50;

do_constraint_plot     = true;
do_cruise_thrust_check = true;
do_ceiling_check       = true;

M_ceiling         = Mach_cruise;
h_ceiling_max_ft  = 65000;
h_ceiling_step_ft = 250;

%% SECTION 2 MISSION WEIGHT CHAIN

fprintf('=== ASW AIRCRAFT CONCEPTUAL SIZING — RAYMER METHOD ===\n\n');

V_cruise_kt = mach_to_tas_knots(Mach_cruise, h_cruise_ft);
t_leg_hr    = Range_leg_nmi / V_cruise_kt;

% Fixed segment fractions (Raymer Table 3.2)
W1_W0 = 0.990; W2_W1 = 0.990; W3_W2 = 0.980;
W7_W6 = 0.990; W8_W7 = 0.992;

W4_W3 = exp(-(t_leg_hr     * c_cruise) / LD_cruise);
W5_W4 = exp(-(E_loiter_hr  * c_loiter) / LD_loiter);
W6_W5 = W4_W3;
W9_W8 = exp(-(E_reserve_hr * c_loiter) / LD_loiter);

Mff_nom   = W1_W0*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5*W7_W6*W8_W7*W9_W8;
Wf_W0_nom = 1 - Mff_nom;
Wf_W0     = Wf_W0_nom + fuel_cont * (1 - Wf_W0_nom);

Wf_W0_floor = 0.42;   % Historical floor: P-3C ~0.44, P-8A ~0.42
if Wf_W0 < Wf_W0_floor
    fprintf('NOTE: Wf/W0 = %.4f below patrol-class floor %.2f — applying floor.\n\n', Wf_W0, Wf_W0_floor);
    Wf_W0 = Wf_W0_floor;
end

WmidCruise_over_W0 = W1_W0 * W2_W1 * W3_W2 * sqrt(W4_W3);

fprintf('--- Mission Segment Fuel Fractions ---\n');
fprintf('  Warmup   W1/W0 : %.4f\n', W1_W0);
fprintf('  Taxi/TO  W2/W1 : %.4f\n', W2_W1);
fprintf('  Climb    W3/W2 : %.4f\n', W3_W2);
fprintf('  Cruise   W4/W3 : %.4f\n', W4_W3);
fprintf('  Loiter   W5/W4 : %.4f\n', W5_W4);
fprintf('  Return   W6/W5 : %.4f  [munitions dropped at W5]\n', W6_W5);
fprintf('  Descent  W7/W6 : %.4f\n', W7_W6);
fprintf('  Land     W8/W7 : %.4f\n', W8_W7);
fprintf('  Reserve  W9/W8 : %.4f\n', W9_W8);
fprintf('  Mff = %.4f   |   Wf/W0 (%.0f%% cont.) = %.4f\n\n', Mff_nom, fuel_cont*100, Wf_W0);

%% SECTION 3 ITERATIVE SIZING LOOP

for outer = 1:100

    W0   = solve_W0(W_payload, Wf_W0, A_ew, C_ew, W0_guess, tol);
    W0_N = lb_to_N(W0);

    sigma_cruise = density_ratio_ISA(h_cruise_ft);
    rho_cruise   = rhoSL * sigma_cruise;
    V_cruise_ms  = V_cruise_kt / kt_per_ms;
    q_cruise     = 0.5 * rho_cruise * V_cruise_ms^2;

    S_ft2 = 1680;   % Fixed wing area (ft^2)
    S_m2  = S_ft2 * 0.09290304;

    WS_MTO_lbft2    = W0 / S_ft2;
    WS_MTO_Nm2      = W0_N / S_m2;
    WS_design_lbft2 = WS_MTO_lbft2;

    TW_req_takeoff = WS_MTO_lbft2 / (CL_TO * TOP_jet);
    TW_design      = thrust_margin * TW_req_takeoff;
    T_total_lbf    = TW_design * W0;
    T_each_lbf     = T_total_lbf / nEng;

    TOFL_actual_ft = TOP_jet * WS_MTO_lbft2 / CLmax_TO;
    takeoff_ok     = TOFL_actual_ft <= TOFL_req_ft;

    rhoSL_slug  = 0.002377;
    g_fps2      = 32.174;
    gamma_land  = 3.0 * pi / 180;
    mu_eff_land = 0.40;

    WS_land_lbft2 = 0;
    for WS_try = 300:-0.5:10
        Vs_try    = sqrt(2 * WS_try * W_land_over_W0 / (rhoSL_slug * CLmax_L));
        Vapp_try  = 1.30 * Vs_try;
        Vtd_try   = 0.95 * Vapp_try;
        s_air_t   = 50 / tan(gamma_land);
        s_flare_t = (Vapp_try^2 / (0.2 * g_fps2)) * sin(gamma_land);
        s_grnd_t  = Vtd_try^2 / (2 * g_fps2 * mu_eff_land);
        if (s_air_t + s_flare_t + s_grnd_t) <= LFL_req_ft
            WS_land_lbft2 = WS_try;
            break;
        end
    end

    Vs_L_land    = sqrt(2 * WS_MTO_lbft2 * W_land_over_W0 / (rhoSL_slug * CLmax_L));
    Vapp_land    = 1.30 * Vs_L_land;
    s_air_act    = 50 / tan(gamma_land);
    s_flare_act  = (Vapp_land^2 / (0.2 * g_fps2)) * sin(gamma_land);
    s_grnd_act   = (0.95*Vapp_land)^2 / (2 * g_fps2 * mu_eff_land);
    LFL_actual_ft = s_air_act + s_flare_act + s_grnd_act;
    Vapp_land_kt  = Vapp_land / 1.68781;
    Vs_L_land_kt  = Vs_L_land / 1.68781;
    landing_ok    = WS_MTO_lbft2 <= WS_land_lbft2;

    WS_cruiseCap_Nm2   = (q_cruise * CL_design_SC) / max(WmidCruise_over_W0, 1e-6);
    WS_cruiseCap_lbft2 = Nm2_to_lbft2(WS_cruiseCap_Nm2);
    CL_cruise_check    = (WS_MTO_Nm2 * WmidCruise_over_W0) / q_cruise;
    cruiseCL_ok        = CL_cruise_check <= (CL_design_SC + 0.10);

    outer_err = abs(W0 - W0_guess) / max(W0_guess, 1);
    W0_guess  = W0;
    if outer_err < 5e-5; break; end
end

%% SECTION 4 FINAL WEIGHTS & SEGMENT CHAIN

We = (A_ew * W0^C_ew) * W0;
Wf = Wf_W0 * W0;

W1 = W0*W1_W0; W2 = W1*W2_W1; W3 = W2*W3_W2;
W4 = W3*W4_W3; W5 = W4*W5_W4;

W5_corrected = W5 - W_mun;   % Munitions dropped after loiter

W6 = W5_corrected*W6_W5;
W7 = W6*W7_W6; W8 = W7*W8_W7; W9 = W8*W9_W8;

Fuel_warmup_lb   = W0 - W1;
Fuel_taxi_lb     = W1 - W2;
Fuel_climb_lb    = W2 - W3;
Fuel_outbound_lb = W3 - W4;
Fuel_loiter_lb   = W4 - W5;
Fuel_return_lb   = W5_corrected - W6;
Fuel_descent_lb  = W6 - W7;
Fuel_land_lb     = W7 - W8;
Fuel_reserve_lb  = W8 - W9;

Fuel_total_check = Fuel_warmup_lb + Fuel_taxi_lb + Fuel_climb_lb + ...
                   Fuel_outbound_lb + Fuel_loiter_lb + Fuel_return_lb + ...
                   Fuel_descent_lb + Fuel_land_lb + Fuel_reserve_lb;

engine_ok_SL = (T_each_lbf <= LEAP_rating_lbf_SL);

%% SECTION 5 CRUISE THRUST CHECK

if do_cruise_thrust_check
    W_cruise_lbf     = WmidCruise_over_W0 * W0;
    T_req_each_lbf   = (W_cruise_lbf / LD_cruise) / nEng;
    lapse_cruise     = thrust_lapse_turbofan(h_cruise_ft, Mach_cruise);
    T_avail_each_lbf = LEAP_rating_lbf_SL * lapse_cruise;
    engine_ok_cruise = (T_avail_each_lbf >= 1.05 * T_req_each_lbf);
else
    lapse_cruise     = NaN;
    T_req_each_lbf   = NaN;
    T_avail_each_lbf = NaN;
    engine_ok_cruise = true;
end

%% SECTION 6 SPEEDS

rhoSL_slug = 0.002377;
Vs_L_fps   = sqrt(2 * WS_MTO_lbft2 * W_land_over_W0 / (rhoSL_slug * CLmax_L));
Vapp_fps   = 1.30 * Vs_L_fps;
Vs_TO_fps  = sqrt(2 * WS_MTO_lbft2 / (rhoSL_slug * CLmax_TO));
V2_fps     = 1.20 * Vs_TO_fps;

Vs_L_kt  = Vs_L_fps  / 1.68781;
Vapp_kt  = Vapp_fps  / 1.68781;
Vs_TO_kt = Vs_TO_fps / 1.68781;
V2_kt    = V2_fps    / 1.68781;
V2_ms    = V2_fps * 0.3048;

%% SECTION 7 — CRUISE AERODYNAMIC CONSISTENCY

W_cruise_N    = lb_to_N(WmidCruise_over_W0 * W0);
WS_cruise_Nm2 = W_cruise_N / S_m2;
CL_cruise_req = WS_cruise_Nm2 / q_cruise;
CD_polar_check = CD0_clean + k_induced * CL_cruise_req^2;
LD_polar_check = CL_cruise_req / CD_polar_check;

TW_cruise_req       = (q_cruise*CD0_clean)/WS_cruise_Nm2 + k_induced*WS_cruise_Nm2/q_cruise;
TW_cruise_req_SL_eq = TW_cruise_req / max(lapse_cruise, 1e-6);

%% SECTION 8 EFFICIENCY METRICS

Range_total_nmi       = 2 * Range_leg_nmi;
specific_range_nmi_lb = Range_total_nmi / max(Wf, 1);
payload_range_eff     = W_payload * Range_total_nmi;
transport_eff         = payload_range_eff / max(Wf, 1);
breguet_SR            = V_cruise_kt * LD_cruise / c_cruise;
endurance_factor      = LD_loiter / c_loiter;

%% SECTION 9 RESULTS PRINTOUT

fprintf('--- Mission Parameters ---\n');
fprintf('  Cruise M%.2f / %.0f ft  |  TAS = %.1f kt  |  Leg = %.2f hr\n', Mach_cruise, h_cruise_ft, V_cruise_kt, t_leg_hr);
fprintf('  Range = %.0f nmi  |  Loiter = %.1f hr  |  Reserve = %.1f hr\n', Range_total_nmi, E_loiter_hr, E_reserve_hr);
fprintf('  TSFC = %.3f / %.3f hr^-1  |  L/D = %.1f / %.1f\n\n', c_cruise, c_loiter, LD_cruise, LD_loiter);

fprintf('--- Weight Summary ---\n');
fprintf('  MTOW  W0  : %.0f lb  (%.2f tons)\n',       W0, W0/2000);
fprintf('  Empty We  : %.0f lb  (%.1f%% of W0)\n',    We, 100*We/W0);
fprintf('  Fuel  Wf  : %.0f lb  (%.1f%%) [fraction]  |  %.0f lb [segment sum]\n', Wf, 100*Wf/W0, Fuel_total_check);
fprintf('  Payload   : %.0f lb  (%.1f%% of W0)\n',    W_payload, 100*W_payload/W0);
fprintf('  Closure   : %.0f lb  (should = W0)\n\n',   We + Wf + W_payload);

fprintf('--- Wing & Field Performance ---\n');
fprintf('  S = %.0f ft^2 (fixed)  |  b = %.1f ft  |  W/S = %.1f lb/ft^2\n', S_ft2, sqrt(AR_wing*S_ft2), WS_MTO_lbft2);
fprintf('  TOFL = %.0f ft  (req <= %.0f) — %s\n', TOFL_actual_ft, TOFL_req_ft, pass_fail(takeoff_ok));
fprintf('  LFL  = %.0f ft  (s_air=%.0f, s_flare=%.0f, s_gr=%.0f)\n', LFL_actual_ft, s_air_act, s_flare_act, s_grnd_act);
fprintf('  Landing W/S limit = %.1f lb/ft^2 — %s\n', WS_land_lbft2, pass_fail(landing_ok));
fprintf('  Vs_L = %.1f kt  |  Vapp = %.1f kt\n', Vs_L_land_kt, Vapp_land_kt);
fprintf('  Cruise CL = %.4f (target %.2f) — %s\n\n', CL_cruise_check, CL_design_SC, pass_fail(cruiseCL_ok));

fprintf('--- Thrust Sizing ---\n');
fprintf('  T/W req = %.4f  |  T/W design = %.4f\n', TW_req_takeoff, TW_design);
fprintf('  Total SL thrust = %.0f lbf  |  Per engine = %.0f lbf\n', T_total_lbf, T_each_lbf);
fprintf('  LEAP-1A26 rating = %.0f lbf  |  Margin = %.1f%% — %s\n\n', LEAP_rating_lbf_SL, 100*(LEAP_rating_lbf_SL-T_each_lbf)/T_each_lbf, pass_fail(engine_ok_SL));

fprintf('--- Cruise Thrust Check ---\n');
fprintf('  Mid-cruise W/W0 = %.4f  |  Lapse = %.4f\n', WmidCruise_over_W0, lapse_cruise);
fprintf('  T_req = %.0f lbf  |  T_avail = %.0f lbf — %s\n\n', T_req_each_lbf, T_avail_each_lbf, pass_fail(engine_ok_cruise));

fprintf('--- Speeds ---\n');
fprintf('  Vs_L = %.1f kt  |  Vapp = %.1f kt  |  Vs_TO = %.1f kt  |  V2 = %.1f kt\n\n', Vs_L_kt, Vapp_kt, Vs_TO_kt, V2_kt);

fprintf('--- Cruise Aero Consistency ---\n');
fprintf('  CL = %.4f  |  CD polar = %.5f  |  L/D polar = %.2f (assumed %.1f)\n', CL_cruise_req, CD_polar_check, LD_polar_check, LD_cruise);
fprintf('  Cruise T/W (SL equiv) = %.4f  |  Iterations = %d\n\n', TW_cruise_req_SL_eq, outer);

fprintf('--- Fuel by Segment ---\n');
fprintf('  Warmup / Taxi     : %.0f / %.0f lb\n', Fuel_warmup_lb, Fuel_taxi_lb);
fprintf('  Climb             : %.0f lb\n', Fuel_climb_lb);
fprintf('  Outbound cruise   : %.0f lb\n', Fuel_outbound_lb);
fprintf('  Loiter            : %.0f lb\n', Fuel_loiter_lb);
fprintf('  [Mun. drop: %.0f lb — not fuel]\n', W_mun);
fprintf('  Return cruise     : %.0f lb\n', Fuel_return_lb);
fprintf('  Descent / Land    : %.0f / %.0f lb\n', Fuel_descent_lb, Fuel_land_lb);
fprintf('  Reserve           : %.0f lb\n', Fuel_reserve_lb);
fprintf('  Total             : %.0f lb\n\n', Fuel_total_check);

fprintf('--- Efficiency ---\n');
fprintf('  Specific range = %.4f nmi/lb  |  Breguet SR = %.1f nmi\n', specific_range_nmi_lb, breguet_SR);
fprintf('  Endurance factor = %.2f hr  |  Transport eff = %.2f payload-nmi/lb\n\n', endurance_factor, transport_eff);

%% SECTION 10 CONSTRAINT DIAGRAM

WS_sweep_lbft2 = linspace(20, 220, 600);
WS_sweep_Nm2   = WS_sweep_lbft2 * 47.88025898;

TW_takeoff = WS_sweep_lbft2 ./ (CL_TO * TOP_jet);

TW_cruise_arr = zeros(size(WS_sweep_Nm2));
for i = 1:length(WS_sweep_Nm2)
    WS_i = WS_sweep_Nm2(i);
    TW_cruise_arr(i) = (q_cruise*CD0_clean)/WS_i + k_induced*WS_i/q_cruise;
end
TW_cruise_SL = TW_cruise_arr / max(lapse_cruise, 1e-6);

ROC_ms    = 1500 * 0.00508;
Vclimb_ms = max(V2_ms, 60);
q_climb   = 0.5 * rhoSL * Vclimb_ms^2;

TW_climb_SL = zeros(size(WS_sweep_Nm2));
for i = 1:length(WS_sweep_Nm2)
    WS_i  = WS_sweep_Nm2(i);
    CL_cl = WS_i / q_climb;
    CD_cl = CD0_clean + k_induced * CL_cl^2;
    TW_climb_SL(i) = CD_cl/CL_cl + ROC_ms/Vclimb_ms;
end

if do_constraint_plot
    figure('Name','Constraint Diagram');
    hold on; grid on; box on;
    plot(WS_sweep_lbft2, TW_takeoff,   'b-', 'LineWidth', 2, 'DisplayName', 'Takeoff (6000 ft)');
    plot(WS_sweep_lbft2, TW_cruise_SL, 'r-', 'LineWidth', 2, 'DisplayName', 'Cruise (SL equiv)');
    plot(WS_sweep_lbft2, TW_climb_SL,  'g-', 'LineWidth', 2, 'DisplayName', 'Climb (1500 fpm)');
    % Vertical limit lines — HandleVisibility off, labeled via text
    xline(WS_land_lbft2,      '-.', 'Color', [0.8 0 0.8],  'LineWidth', 1.8, 'HandleVisibility', 'off');
    xline(WS_cruiseCap_lbft2, '--', 'Color', [0 0.75 0.75], 'LineWidth', 1.8, 'HandleVisibility', 'off');
    xline(WS_MTO_lbft2,       '--k', 'LineWidth', 1.5,      'HandleVisibility', 'off');
    text(WS_land_lbft2+0.8,      0.04, 'Landing W/S limit',       'Color', [0.8 0 0.8],   'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'bottom');
    text(WS_cruiseCap_lbft2+0.8, 0.04, 'Cruise C_L cap',          'Color', [0 0.75 0.75], 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'bottom');
    text(WS_MTO_lbft2+0.8,       0.04, sprintf('Design W/S = %.0f lb/ft^2', WS_MTO_lbft2), ...
         'Color', 'k', 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'bottom');
    % Dummy handles for legend
    plot(NaN, NaN, '-.', 'Color', [0.8 0 0.8],   'LineWidth', 1.5, 'DisplayName', sprintf('Landing W/S limit  %.0f lb/ft^2', WS_land_lbft2));
    plot(NaN, NaN, '--', 'Color', [0 0.75 0.75],  'LineWidth', 1.5, 'DisplayName', sprintf('Cruise C_L cap  %.0f lb/ft^2', WS_cruiseCap_lbft2));
    plot(NaN, NaN, '--k', 'LineWidth', 1.5,        'DisplayName', sprintf('Design W/S  %.0f lb/ft^2', WS_MTO_lbft2));
    plot(WS_MTO_lbft2, TW_design, 'ko', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'DisplayName', 'Design Point');
    xlabel('Wing Loading  W/S  (lb/ft^2)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Thrust-to-Weight  T/W  (SL)',  'FontSize', 12, 'FontWeight', 'bold');
    title('ASW Aircraft Constraint Diagram', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10);
    ylim([0, 0.7]); xlim([20, 220]);
end

%% SECTION 11 CEILING ANALYSIS

if do_ceiling_check

    W_ceiling_lbf = WmidCruise_over_W0 * W0;   % Conservative: outbound mid-cruise
    W_ceiling_N   = lb_to_N(W_ceiling_lbf);
    CLmax_cruise  = 0.87;   % Buffet onset limit at M0.75

    h_grid_ft = 0 : h_ceiling_step_ft : h_ceiling_max_ft;
    T_avail   = zeros(size(h_grid_ft));
    T_req     = zeros(size(h_grid_ft));
    ROC_fpm   = zeros(size(h_grid_ft));

    for i = 1:length(h_grid_ft)
        h_ft  = h_grid_ft(i);
        rho_i = rhoSL * density_ratio_ISA(h_ft);
        V_ms  = mach_to_tas_knots(M_ceiling, h_ft) / kt_per_ms;
        q_i   = 0.5 * rho_i * V_ms^2;
        CL_i  = W_ceiling_N / (q_i * S_m2);

        if CL_i > CLmax_cruise
            T_req(i) = Inf; T_avail(i) = 0; ROC_fpm(i) = -9999;
            continue
        end

        CD_i       = CD0_clean + k_induced * CL_i^2;
        T_req(i)   = W_ceiling_lbf * CD_i / CL_i;
        T_avail(i) = nEng * LEAP_rating_lbf_SL * thrust_lapse_turbofan(h_ft, M_ceiling);
        ROC_fpm(i) = (T_avail(i) - T_req(i)) / W_ceiling_lbf * V_ms * 196.8504;
    end

    ROC_crit = 100;
    idx = find(ROC_fpm(1:end-1) >= ROC_crit & ROC_fpm(2:end) <= ROC_crit, 1);
    if ~isempty(idx)
        h1 = h_grid_ft(idx); h2 = h_grid_ft(idx+1);
        r1 = ROC_fpm(idx);   r2 = ROC_fpm(idx+1);
        if r2 <= -9000
            service_ceiling_ft = h1;
            ceiling_limit_type = 'aerodynamic (buffet onset)';
        else
            service_ceiling_ft = h1 + (ROC_crit-r1)*(h2-h1)/(r2-r1);
            ceiling_limit_type = 'thrust limited';
        end
    elseif ROC_fpm(end) > ROC_crit
        service_ceiling_ft = h_ceiling_max_ft;
        ceiling_limit_type = 'exceeds scan range';
        warning('Ceiling exceeds scan range of %.0f ft', h_ceiling_max_ft);
    else
        service_ceiling_ft = NaN;
        ceiling_limit_type = 'not found';
    end

    ceiling_ok = ~isnan(service_ceiling_ft) && service_ceiling_ft >= 30000;

    fprintf('--- Ceiling Analysis ---\n');
    fprintf('  Analysis weight (mid-cruise) : %.0f lb\n', W_ceiling_lbf);
    if ~isnan(service_ceiling_ft)
        fprintf('  Service ceiling (100 fpm)    : %.0f ft  [%s]\n', service_ceiling_ft, ceiling_limit_type);
    else
        fprintf('  Service ceiling              : NOT FOUND\n');
    end
    fprintf('  30,000 ft requirement        : %s\n\n', pass_fail(ceiling_ok));

    figure('Name','Thrust vs Altitude');
    hold on; grid on; box on;
    plot(h_grid_ft/1000, T_avail, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Thrust Available');
    plot(h_grid_ft/1000, T_req,   'r-', 'LineWidth', 2.5, 'DisplayName', 'Thrust Required');
    % Reference lines — use text annotations so xline doesn't pollute legend
    xline(30, '-.', 'Color', [0 0.6 0], 'LineWidth', 1.8, 'HandleVisibility', 'off');
    text(30.6, max(T_avail(isfinite(T_avail)))*0.92, '30 kft req', ...
         'Color', [0 0.6 0], 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'top');
    if ~isnan(service_ceiling_ft)
        xline(service_ceiling_ft/1000, '--k', 'LineWidth', 1.8, 'HandleVisibility', 'off');
        text(service_ceiling_ft/1000 + 0.6, max(T_avail(isfinite(T_avail)))*0.92, ...
             sprintf('Ceiling %.0f ft', service_ceiling_ft), ...
             'Color', 'k', 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'top');
        plot(NaN, NaN, '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Service ceiling %.0f ft', service_ceiling_ft));
    end
    plot(NaN, NaN, '-.', 'Color', [0 0.6 0], 'LineWidth', 1.5, 'DisplayName', '30,000 ft requirement');
    xlabel('Altitude (1000 ft)','FontSize',12,'FontWeight','bold');
    ylabel('Thrust (lbf)',      'FontSize',12,'FontWeight','bold');
    title('Available vs Required Thrust — Ceiling Analysis','FontSize',13,'FontWeight','bold');
    legend('Location','northeast','FontSize',10);

    % Clip ROC so the aerodynamic cutoff shows as a clean drop to NaN, not -9999
    ROC_plot = ROC_fpm;
    ROC_plot(ROC_plot <= -9000) = NaN;

    figure('Name','Rate of Climb vs Altitude');
    hold on; grid on; box on;
    plot(h_grid_ft/1000, ROC_plot, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Rate of Climb (mid-cruise weight)');
    % Reference lines with HandleVisibility off — use dummy plots for legend
    yline(ROC_crit, '--r', 'LineWidth', 1.8, 'HandleVisibility', 'off');
    xline(30, '-.', 'Color', [0 0.6 0], 'LineWidth', 1.8, 'HandleVisibility', 'off');
    text(30.6, 9200, '30 kft req', 'Color', [0 0.6 0], 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'top');
    if ~isnan(service_ceiling_ft)
        xline(service_ceiling_ft/1000, '--k', 'LineWidth', 1.8, 'HandleVisibility', 'off');
        text(service_ceiling_ft/1000 + 0.6, 9200, sprintf('Ceiling %.0f ft', service_ceiling_ft), ...
             'Color', 'k', 'FontSize', 9, 'Rotation', 90, 'VerticalAlignment', 'top');
        plot(NaN, NaN, '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Service ceiling %.0f ft', service_ceiling_ft));
    end
    plot(NaN, NaN, '--r',  'LineWidth', 1.5, 'DisplayName', '100 fpm — service ceiling threshold');
    plot(NaN, NaN, '-.', 'Color', [0 0.6 0], 'LineWidth', 1.5, 'DisplayName', '30,000 ft requirement');
    xlabel('Altitude (1000 ft)', 'FontSize',12,'FontWeight','bold');
    ylabel('Rate of Climb (fpm)','FontSize',12,'FontWeight','bold');
    title('Rate of Climb vs Altitude — Ceiling Analysis','FontSize',13,'FontWeight','bold');
    legend('Location','northeast','FontSize',10);
    ylim([-500, 10000]);
end

%% SECTION 12 WEIGHT BREAKDOWN PIE CHART

figure('Name','Weight Breakdown');
pie([We, Wf, W_payload]);
legend({sprintf('Empty  %.0f lb (%.0f%%)',   We,        round(100*We/W0)), ...
        sprintf('Fuel   %.0f lb (%.0f%%)',    Wf,        round(100*Wf/W0)), ...
        sprintf('Payload %.0f lb (%.0f%%)',   W_payload, round(100*W_payload/W0))}, ...
       'Location','southoutside','Orientation','horizontal','FontSize',10);
title(sprintf('Weight Breakdown — MTOW = %.0f lb  (%.2f tons)', W0, W0/2000), ...
      'FontSize',12,'FontWeight','bold');
colormap([0.30 0.50 0.85; 0.90 0.50 0.20; 0.35 0.72 0.40]);

fprintf('=== Sizing complete. All figures docked. ===\n');

%% SECTION 13 RATE OF CLIMB vs ALTITUDE (Multi-Weight)

fprintf('\n--- Rate of Climb Analysis ---\n');

h_roc_ft    = 0 : 500 : h_ceiling_max_ft;
W_cases_lb  = [W0,  WmidCruise_over_W0*W0,  W5_corrected];
W_case_lbls = {'MTOW', 'Mid-cruise (outbound)', 'Post-loiter (mun. dropped)'};
W_case_clrs = {'b', 'r', 'g'};

figure('Name','Vertical Rate of Climb vs Altitude');
hold on; grid on; box on;

for wc = 1:3
    W_case_lbf = W_cases_lb(wc);
    W_case_N   = lb_to_N(W_case_lbf);
    ROC_case   = zeros(size(h_roc_ft));

    for i = 1:length(h_roc_ft)
        h_ft  = h_roc_ft(i);
        rho_i = rhoSL * density_ratio_ISA(h_ft);
        V_ms  = mach_to_tas_knots(Mach_cruise, h_ft) / kt_per_ms;
        q_i   = 0.5 * rho_i * V_ms^2;
        CL_i  = W_case_N / (q_i * S_m2);
        if CL_i > CLmax_cruise; ROC_case(i) = NaN; continue; end
        CD_i        = CD0_clean + k_induced * CL_i^2;
        T_avail_i   = nEng * LEAP_rating_lbf_SL * thrust_lapse_turbofan(h_ft, Mach_cruise);
        ROC_case(i) = (T_avail_i - W_case_lbf*CD_i/CL_i) / W_case_lbf * V_ms * 196.8504;
    end

    valid = ~isnan(ROC_case);
    h_v = h_roc_ft(valid); r_v = ROC_case(valid);
    cidx = find(r_v(1:end-1) >= 100 & r_v(2:end) <= 100, 1);
    if ~isempty(cidx)
        ceil_h = h_v(cidx) + (100-r_v(cidx))*(h_v(cidx+1)-h_v(cidx))/(r_v(cidx+1)-r_v(cidx));
        fprintf('  %-30s: ceiling = %.0f ft  |  SL ROC = %.0f fpm\n', W_case_lbls{wc}, ceil_h, max(ROC_case(1),0));
    else
        fprintf('  %-30s: ceiling > %.0f ft  |  SL ROC = %.0f fpm\n', W_case_lbls{wc}, h_ceiling_max_ft, max(ROC_case(1),0));
    end

    plot(h_roc_ft/1000, ROC_case, [W_case_clrs{wc} '-'], 'LineWidth', 2, 'DisplayName', W_case_lbls{wc});
end

% Horizontal reference lines
yline(100, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8, 'HandleVisibility', 'off');
yline(500, '-.', 'Color', [0.8 0 0.8],   'LineWidth', 1.8, 'HandleVisibility', 'off');
% Cruise altitude vertical line
xline(h_cruise_ft/1000, '--c', 'LineWidth', 1.8, 'HandleVisibility', 'off');
text(h_cruise_ft/1000 + 0.6, 9500, sprintf('Cruise %d kft', h_cruise_ft/1000), ...
     'Color', [0 0.8 0.8], 'FontSize', 10, 'Rotation', 90, 'VerticalAlignment', 'top', 'FontWeight', 'bold');
% Dummy plot handles for legend
plot(NaN, NaN, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5, 'DisplayName', '100 fpm — Service ceiling');
plot(NaN, NaN, '-.', 'Color', [0.8 0 0.8],   'LineWidth', 1.5, 'DisplayName', '500 fpm — Operational minimum');
plot(NaN, NaN, '--c', 'LineWidth', 1.5, 'DisplayName', sprintf('Cruise altitude %d kft', h_cruise_ft/1000));
xlabel('Altitude (1000 ft)', 'FontSize',13,'FontWeight','bold');
ylabel('Rate of Climb (fpm)','FontSize',13,'FontWeight','bold');
title('Rate of Climb vs Altitude — Multi-Weight','FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',10);
xlim([0, h_ceiling_max_ft/1000]); ylim([-500, 10500]);

%% SECTION 14 PAYLOAD-RANGE-ENDURANCE DIAGRAM

fprintf('\n--- Payload-Range-Endurance ---\n');

npts = 500;
f_overhead_frac = 1 - (W1_W0*W2_W1*W3_W2*W7_W6*W8_W7*W9_W8);
f_usable        = 1 - f_overhead_frac;
f_cr = 33289 / (33289 + 12972);
f_lo = 12972 / (33289 + 12972);

Wf_s1 = linspace(500, Wf, npts);
R_s1 = zeros(1,npts); E_s1 = zeros(1,npts);
for k = 1:npts
    W0k     = We + W_payload + Wf_s1(k);
    Wfu_k   = Wf_s1(k) * f_usable;
    Wi_cr   = W0k * W1_W0 * W2_W1 * W3_W2;
    Wf_leg1 = max(Wi_cr - Wfu_k*f_cr/2, Wi_cr*0.60);
    if Wi_cr > Wf_leg1; R_s1(k) = 2*(V_cruise_kt/c_cruise)*LD_cruise*log(Wi_cr/Wf_leg1); end
    Wi_lo = max(Wi_cr - Wfu_k*f_cr, Wi_cr*0.55);
    Wf_lo = max(Wi_lo - Wfu_k*f_lo, Wi_lo*0.75);
    if Wi_lo > Wf_lo; E_s1(k) = (LD_loiter/c_loiter)*log(Wi_lo/Wf_lo); end
end

Wp_s2 = linspace(W_payload, 0, npts);
Wf_s2 = W0 - We - Wp_s2;
R_s2 = zeros(1,npts); E_s2 = zeros(1,npts);
for k = 1:npts
    Wfu_k   = Wf_s2(k) * f_usable;
    Wi_cr   = W0 * W1_W0 * W2_W1 * W3_W2;
    Wf_leg1 = max(Wi_cr - Wfu_k*f_cr/2, Wi_cr*0.55);
    if Wi_cr > Wf_leg1; R_s2(k) = 2*(V_cruise_kt/c_cruise)*LD_cruise*log(Wi_cr/Wf_leg1); end
    Wi_lo = max(Wi_cr - Wfu_k*f_cr, Wi_cr*0.50);
    Wf_lo = max(Wi_lo - Wfu_k*f_lo, Wi_lo*0.70);
    if Wi_lo > Wf_lo; E_s2(k) = (LD_loiter/c_loiter)*log(Wi_lo/Wf_lo); end
end

Wi_ferry = W0 * W1_W0 * W2_W1 * W3_W2;
Wf_end_f = max(Wi_ferry - (W0-We)*f_usable*0.85, Wi_ferry*0.50);
R_ferry  = (V_cruise_kt/c_cruise)*LD_cruise*log(Wi_ferry/Wf_end_f);
E_ferry  = R_ferry / V_cruise_kt;

fprintf('  Corner point: %.0f nmi  |  %.1f hr\n', R_s1(end), E_s1(end));
fprintf('  Ferry range:  %.0f nmi\n\n', R_ferry);

figure('Name','Payload-Range-Endurance');
subplot(1,2,1); hold on; grid on; box on;
plot(R_s1, W_payload*ones(1,npts), 'b-', 'LineWidth', 2.5, 'DisplayName', sprintf('Seg 1: full payload (%.0f lb)', W_payload));
plot(R_s2, Wp_s2, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Seg 2: payload-fuel trade');
plot(R_s1(end), W_payload, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'DisplayName', sprintf('Corner point (%.0f nmi)', R_s1(end)));
xline(2*Range_leg_nmi, '--', 'Color', [0.7 0 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(2*Range_leg_nmi + 80, W_payload*0.55, sprintf('Design\n%.0f nmi', 2*Range_leg_nmi), ...
     'Color', [0.7 0 0.7], 'FontSize', 9, 'FontWeight', 'bold');
plot(NaN, NaN, '--', 'Color', [0.7 0 0.7], 'LineWidth', 1.5, 'DisplayName', sprintf('Design mission %.0f nmi', 2*Range_leg_nmi));
text(200, W_payload*0.08, sprintf('\\rightarrow Ferry range: %.0f nmi', R_ferry), 'FontSize', 9, 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold');
xlabel('Range (nmi)','FontSize',12,'FontWeight','bold');
ylabel('Payload (lb)','FontSize',12,'FontWeight','bold');
title('Payload vs Range','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',9); xlim([0,8000]); ylim([0,W_payload*1.25]);

subplot(1,2,2); hold on; grid on; box on;
plot(E_s1, W_payload*ones(1,npts), 'b-', 'LineWidth', 2.5, 'DisplayName', sprintf('Seg 1: full payload (%.0f lb)', W_payload));
plot(E_s2, Wp_s2, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Seg 2: payload-fuel trade');
plot(E_s1(end), W_payload, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'DisplayName', sprintf('Corner point (%.1f hr)', E_s1(end)));
xline(E_loiter_hr, '--', 'Color', [0.7 0 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(E_loiter_hr + 0.15, W_payload*0.55, sprintf('Design\n%.1f hr', E_loiter_hr), ...
     'Color', [0.7 0 0.7], 'FontSize', 9, 'FontWeight', 'bold');
plot(NaN, NaN, '--', 'Color', [0.7 0 0.7], 'LineWidth', 1.5, 'DisplayName', sprintf('Design loiter %.1f hr', E_loiter_hr));
text(0.2, W_payload*0.08, sprintf('\\rightarrow Ferry endurance: %.1f hr', E_ferry), 'FontSize', 9, 'Color', [0.1 0.6 0.1], 'FontWeight', 'bold');
xlabel('Endurance (hr)','FontSize',12,'FontWeight','bold');
ylabel('Payload (lb)','FontSize',12,'FontWeight','bold');
title('Payload vs Endurance','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',9); xlim([0,16]); ylim([0,W_payload*1.25]);
sgtitle('Payload-Range-Endurance — ASW Aircraft','FontSize',14,'FontWeight','bold');

%% SECTION 15 ENGINE PERFORMANCE vs ALTITUDE

fprintf('--- Engine Performance vs Altitude ---\n');

h_eng_ft    = 0 : 500 : h_ceiling_max_ft;
T_avail_eng = zeros(size(h_eng_ft));
lapse_eng   = zeros(size(h_eng_ft));
TSFC_eng    = zeros(size(h_eng_ft));
SEP_eng     = zeros(size(h_eng_ft));

W_eng_lbf = WmidCruise_over_W0 * W0;
W_eng_N   = lb_to_N(W_eng_lbf);

for i = 1:length(h_eng_ft)
    h_ft          = h_eng_ft(i);
    lapse_i       = thrust_lapse_turbofan(h_ft, Mach_cruise);
    lapse_eng(i)  = lapse_i;
    T_avail_eng(i)= nEng * LEAP_rating_lbf_SL * lapse_i;

    sigma_h     = density_ratio_ISA(h_ft);
    TSFC_eng(i) = c_cruise * (1.0 + 0.15*(1-sigma_h^0.3) + 0.10*max(0,(h_ft-40000)/20000));
    TSFC_eng(i) = max(TSFC_eng(i), c_cruise*0.95);

    V_ms  = mach_to_tas_knots(Mach_cruise, h_ft) / kt_per_ms;
    q_i   = 0.5 * rhoSL * sigma_h * V_ms^2;
    CL_i  = W_eng_N / (q_i * S_m2);
    if CL_i > CLmax_cruise
        SEP_eng(i) = NaN;
    else
        CD_i       = CD0_clean + k_induced * CL_i^2;
        SEP_eng(i) = (T_avail_eng(i) - W_eng_lbf*CD_i/CL_i) / W_eng_lbf * V_ms * 196.8504;
    end
end

alts_print = [0, 10000, 20000, 35000, 40000, 45000];
fprintf('  %-10s  %-14s  %-8s  %-12s  %-10s\n', 'Alt (ft)', 'T_avail (lbf)', 'Lapse', 'TSFC (hr^-1)', 'SEP (fpm)');
for h_p = alts_print
    [~, ip] = min(abs(h_eng_ft - h_p));
    fprintf('  %-10d  %-14.0f  %-8.4f  %-12.4f  %-10.0f\n', h_p, T_avail_eng(ip), lapse_eng(ip), TSFC_eng(ip), SEP_eng(ip));
end

figure('Name','Engine Performance vs Altitude');

subplot(2,2,1); hold on; grid on; box on;
plot(h_eng_ft/1000, T_avail_eng, 'b-', 'LineWidth', 2);
xline(h_cruise_ft/1000,    '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
if ~isnan(service_ceiling_ft)
    xline(service_ceiling_ft/1000, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(service_ceiling_ft/1000 + 0.5, max(T_avail_eng)*0.92, 'Ceiling', ...
         'Color','k','FontSize',8,'Rotation',90,'VerticalAlignment','top');
end
text(h_cruise_ft/1000 + 0.5, max(T_avail_eng)*0.92, 'Cruise alt', ...
     'Color','r','FontSize',8,'Rotation',90,'VerticalAlignment','top');
xlabel('Altitude (1000 ft)','FontSize',10,'FontWeight','bold');
ylabel('Total Thrust (lbf)','FontSize',10,'FontWeight','bold');
title('Thrust Available — Both Engines','FontSize',11,'FontWeight','bold');
xlim([0, h_ceiling_max_ft/1000]);

subplot(2,2,2); hold on; grid on; box on;
plot(h_eng_ft/1000, lapse_eng, 'r-', 'LineWidth', 2);
xline(h_cruise_ft/1000, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(h_cruise_ft/1000 + 0.5, 0.95, 'Cruise alt', ...
     'Color','b','FontSize',8,'Rotation',90,'VerticalAlignment','top');
xlabel('Altitude (1000 ft)','FontSize',10,'FontWeight','bold');
ylabel('Thrust Lapse Ratio','FontSize',10,'FontWeight','bold');
title('Thrust Lapse — LEAP-1A26 (BPR=11)','FontSize',11,'FontWeight','bold');
xlim([0, h_ceiling_max_ft/1000]); ylim([0, 1.05]);

subplot(2,2,3); hold on; grid on; box on;
plot(h_eng_ft/1000, TSFC_eng, 'g-', 'LineWidth', 2);
yline(c_cruise, '--b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(h_cruise_ft/1000, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(h_cruise_ft/1000 + 0.5, c_cruise*1.32, 'Cruise alt', ...
     'Color','r','FontSize',8,'Rotation',90,'VerticalAlignment','top');
text(h_ceiling_max_ft/1000*0.55, c_cruise*1.005, sprintf('Design TSFC = %.2f hr^{-1}', c_cruise), ...
     'Color','b','FontSize',8);
xlabel('Altitude (1000 ft)','FontSize',10,'FontWeight','bold');
ylabel('TSFC (hr^{-1})','FontSize',10,'FontWeight','bold');
title('TSFC Estimate vs Altitude','FontSize',11,'FontWeight','bold');
xlim([0, h_ceiling_max_ft/1000]); ylim([c_cruise*0.93, c_cruise*1.35]);

subplot(2,2,4); hold on; grid on; box on;
SEP_valid_max = max(SEP_eng(~isnan(SEP_eng)));
plot(h_eng_ft/1000, SEP_eng, 'm-', 'LineWidth', 2, 'DisplayName', 'SEP (mid-cruise weight)');
% 100 fpm line — use dark red, not yellow (invisible on white)
yline(100, '--', 'Color', [0.7 0 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(h_cruise_ft/1000, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
if ~isnan(service_ceiling_ft)
    xline(service_ceiling_ft/1000, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(service_ceiling_ft/1000 + 0.5, SEP_valid_max*0.90, sprintf('Ceiling %.0f ft', service_ceiling_ft), ...
         'Color','k','FontSize',8,'Rotation',90,'VerticalAlignment','top');
end
text(h_cruise_ft/1000 + 0.5, SEP_valid_max*0.90, 'Cruise alt', ...
     'Color','r','FontSize',8,'Rotation',90,'VerticalAlignment','top');
text(h_ceiling_max_ft/1000*0.35, 200, '100 fpm — service ceiling', ...
     'Color',[0.7 0 0],'FontSize',8);
% Dummy handle
plot(NaN, NaN, '--', 'Color', [0.7 0 0], 'LineWidth', 1.5, 'DisplayName', '100 fpm — service ceiling');
xlabel('Altitude (1000 ft)','FontSize',10,'FontWeight','bold');
ylabel('Specific Excess Power (fpm)','FontSize',10,'FontWeight','bold');
title('SEP at Mid-Cruise Weight','FontSize',11,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
xlim([0, h_ceiling_max_ft/1000]); ylim([-200, SEP_valid_max*1.10]);

sgtitle('CFM LEAP-1A26 Engine Performance vs Altitude (M0.75)','FontSize',13,'FontWeight','bold');

%% HELPER FUNCTIONS

function W0 = solve_W0(W_fixed, Wf_W0, A, C, W0_init, tol)
    W0 = W0_init;
    for it = 1:2000
        We_W0 = A * W0^C;
        denom = 1 - We_W0 - Wf_W0;
        if denom <= 0.02; error('Weight fractions infeasible: denom = %.4f', denom); end
        W0_new = W_fixed / denom;
        if abs(W0_new - W0) / max(W0, 1) < tol; W0 = W0_new; break; end
        W0 = W0_new;
    end
end

function V = mach_to_tas_knots(M, h_ft)
    h = h_ft * 0.3048;
    if h <= 11000; T = 288.15 - 0.0065*h; else; T = 288.15 - 0.0065*11000; end
    V = M * sqrt(1.4 * 287.05 * T) * 1.94384;
end

function lapse = thrust_lapse_turbofan(h_ft, M)
    BPR = 11.0;
    n   = 0.6  + 0.012 * BPR;   % 0.732
    m   = 0.30 - 0.016 * BPR;   % 0.124
    lapse = max(min((density_ratio_ISA(h_ft)^n) * (1 - m*M), 1.0), 0.02);
end

function sigma = density_ratio_ISA(h_ft)
    T0=288.15; p0=101325; g0=9.80665; R=287.058;
    h=h_ft*0.3048;
    h1=11000; h2=20000; h3=32000;
    aL1=-0.0065; aL3=0.001;
    T1=T0+aL1*h1; p1=p0*(T1/T0)^(-g0/(aL1*R));
    T2=T1;        p2=p1*exp(-g0*(h2-h1)/(R*T1));
    T3=T2+aL3*(h3-h2); p3=p2*(T3/T2)^(-g0/(aL3*R));
    if     h<=h1; T=T0+aL1*h;       p=p0*(T/T0)^(-g0/(aL1*R));
    elseif h<=h2; T=T1;              p=p1*exp(-g0*(h-h1)/(R*T1));
    elseif h<=h3; T=T2+aL3*(h-h2);  p=p2*(T/T2)^(-g0/(aL3*R));
    else;         T=T3;              p=p3*exp(-g0*(h-h3)/(R*T3));
    end
    sigma = (p/(R*T)) / (p0/(R*T0));
end

function N     = lb_to_N(lb);      N     = lb * 4.4482216152605; end
function lbft2 = Nm2_to_lbft2(Nm2); lbft2 = Nm2 / 47.88025898;  end
function str   = pass_fail(cond);  if cond; str='PASS'; else; str='*** FAIL ***'; end; end