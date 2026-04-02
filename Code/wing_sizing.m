function out = wing_sizing()
% based on Raymer
% Marzuk Hasan
clc; close all;

%% VARIABLES
% wt.W0_lbf      takeoff gross weight [lbf]
% wt.We_lbf      empty weight [lbf]
% wt.Wf_lbf      fuel weight [lbf]
% wt.Wp_lbf      payload weight [lbf]
% wt.Wc_lbf      crew weight [lbf]
%
% mis.rng_nmi    one-way mission range [nmi]
% mis.lo_hr      loiter time [hr]
% mis.re_hr      reserve loiter time [hr]
% mis.cont       fuel contingency fraction [-]
% mis.LD_cr      cruise lift-to-drag ratio [-]
% mis.LD_lo      loiter lift-to-drag ratio [-]
% mis.c_cr_hr    cruise TSFC [1/hr]
% mis.c_lo_hr    loiter TSFC [1/hr]
%
% cr.M           cruise Mach number [-]
% cr.h_ft        cruise altitude [ft]
% cr.CLdes       design cruise CL target [-]
%
% to.TOFL_ft     takeoff field length requirement [ft]
% ld.LFL_ft      landing field length requirement [ft]
% ld.mL_mTO      landing-to-takeoff mass fraction [-]
% ld.CLmax       landing max lift coefficient [-]
%
% to.CLmax       takeoff max lift coefficient [-]
% to.CLg_frac    takeoff ground-roll lift fraction of CLmax [-]
% to.CD0         takeoff-config parasite drag coefficient [-]
% to.e           takeoff Oswald efficiency [-]
%
% wg.AR          wing aspect ratio [-]
% wg.e           cruise Oswald efficiency [-]
% wg.CD0cr       cruise parasite drag coefficient [-]
%
% eng.n          engine count [-]
% eng.Tsl_lbf    sea-level static rating per engine [lbf]
% eng.mar        thrust margin factor [-]
%
% aero.A         empty-weight regression coefficient [-]
% aero.C         empty-weight regression exponent [-]
%
% atm.rho0       sea-level density [slug/ft^3]
% atm.mu0        sea-level dynamic viscosity [slug/(ft*s)]
% atm.g0         gravity [ft/s^2]
%
% ws.WS_ref      final wing loading [lbf/ft^2]
% ws.Sref        final wing area [ft^2]
% ws.b           wing span [ft]
% ws.cbar        mean aerodynamic chord [ft]
%
% to.TW_avail    installed sea-level static T/W [-]
% to.TW_allow    design allowed T/W after margin [-]
% to.TW_req      takeoff-required T/W [-]
%
% ld.WS_max      landing wing-loading cap [lbf/ft^2]
% to.WS_max      takeoff wing-loading cap [lbf/ft^2]
% cr.WS_max      cruise wing-loading cap [lbf/ft^2]

%% INPUTS
cfg = struct();

cfg.wt.Wp_lbf = 10000 + 20000;
cfg.wt.Wc_lbf = 800;

cfg.mis.rng_nmi = 2000;
cfg.mis.lo_hr    = 4.0;
cfg.mis.re_hr    = 0.5;
cfg.mis.cont     = 0.05;
cfg.mis.LD_cr    = 20;
cfg.mis.LD_lo    = 20;
cfg.mis.c_cr_hr  = 0.55;
cfg.mis.c_lo_hr  = 0.50;

cfg.cr.M      = 0.75;
cfg.cr.h_ft   = 35000;
cfg.cr.CLdes  = 0.40;

cfg.to.TOFL_ft   = 6000;
cfg.ld.LFL_ft    = 5000;
cfg.ld.mL_mTO    = 0.78;
cfg.ld.CLmax     = 2.0;

cfg.to.CLmax     = 2.2;
cfg.to.CLg_frac   = 0.80;   % starting point: CLg = 0.8*CLmax_TO
cfg.to.CD0        = 0.030;   % takeoff-config CD0; tune this
cfg.to.e          = 0.80;    % takeoff-config e; tune this

cfg.wg.AR     = 11.5;
cfg.wg.e      = 0.78;
cfg.wg.CD0cr  = 0.022;

cfg.eng.n       = 2;
cfg.eng.Tsl_lbf = 30000;
cfg.eng.mar     = 1.20;

cfg.aero.A = 0.93;
cfg.aero.C = -0.07;

cfg.atm.g0  = 32.174;


cfg.opt.W0_guess_lbf = 170000;
cfg.opt.tol          = 1e-4;
cfg.opt.doPlot       = true;

%% MISSION FUEL FRACTION
cru.V_fts = isa_speed_of_sound_ftps(cfg.cr.h_ft) * cfg.cr.M;
cru.V_kt  = cru.V_fts / 1.68780986;

cfg.wt.Wfix_lbf = cfg.wt.Wp_lbf + cfg.wt.Wc_lbf;

for it = 1:40
    t_leg_hr = cfg.mis.rng_nmi / cru.V_kt;

    w1_w0 = 0.970;
    w2_w1 = 0.985;
    w3_w2 = 0.995;
    w7_w6 = 0.995;
    w8_w7 = 0.980;

    w4_w3 = exp(-(t_leg_hr * cfg.mis.c_cr_hr) / cfg.mis.LD_cr);
    w5_w4 = exp(-(cfg.mis.lo_hr * cfg.mis.c_lo_hr) / cfg.mis.LD_lo);
    w6_w5 = exp(-(t_leg_hr * cfg.mis.c_cr_hr) / cfg.mis.LD_cr);
    w9_w8 = exp(-(cfg.mis.re_hr * cfg.mis.c_lo_hr) / cfg.mis.LD_lo);

    Mff_nom = w1_w0 * w2_w1 * w3_w2 * w4_w3 * w5_w4 * w6_w5 * w7_w6 * w8_w7 * w9_w8;
    Wf_W0_nom = 1 - Mff_nom;
    Wf_W0     = 1 - (1 - Wf_W0_nom) * (1 - cfg.mis.cont);

    W0_lbf = solve_W0_lbf(cfg.wt.Wfix_lbf, Wf_W0, cfg.aero.A, cfg.aero.C, ...
                          cfg.opt.W0_guess_lbf, cfg.opt.tol);

    err = abs(W0_lbf - cfg.opt.W0_guess_lbf) / max(cfg.opt.W0_guess_lbf, 1);
    cfg.opt.W0_guess_lbf = W0_lbf;

    if err < 5e-4
        break;
    end
end

%% FINAL WEIGHTS
res.wt.W0_lbf = W0_lbf;
res.wt.We_lbf = (cfg.aero.A * W0_lbf^cfg.aero.C) * W0_lbf;
res.wt.Wf_lbf = Wf_W0 * W0_lbf;
res.wt.Wp_lbf = cfg.wt.Wp_lbf;
res.wt.Wc_lbf = cfg.wt.Wc_lbf;
res.wt.Wfix_lbf = cfg.wt.Wfix_lbf;
res.wt.WmidCruise_W0 = w1_w0 * w2_w1 * w3_w2 * sqrt(w4_w3);

%% ATMOSPHERE
% Uses your calcISA(h_ft, Lref_ft, V_fts) function for density, viscosity, Re.
[to.rho_slugft3, to.mu_slugfts, ~] = calcISA(0);
[ld.rho_slugft3, ld.mu_slugfts, ~] = calcISA(0);
[cr.rho_slugft3, cr.mu_slugfts, ~] = calcISA(cfg.cr.h_ft);

cr.a_fts = isa_speed_of_sound_ftps(cfg.cr.h_ft);
cr.V_fts = cr.a_fts * cfg.cr.M;
cr.V_kt   = cr.V_fts / 1.68780986;

to.rho_slugft3 = to.rho_slugft3;
ld.rho_slugft3 = ld.rho_slugft3;

%% BASIC AIRFRAME FACTORS
to.CLg = cfg.to.CLg_frac * cfg.to.CLmax;
to.K   = 1 / (pi * cfg.to.e * cfg.wg.AR);
cr.K   = 1 / (pi * cfg.wg.e * cfg.wg.AR);

res.eng.TW_avail_SL = (cfg.eng.n * cfg.eng.Tsl_lbf) / res.wt.W0_lbf;
res.eng.TW_allow_SL = res.eng.TW_avail_SL / cfg.eng.mar;

%% WING-LOADING CAPS
% Takeoff cap: Raymer ground-roll equation solved for max W/S given allowed T/W
to.WS_max_lbf_ft2 = WS_from_takeoff_Raymer_lbft2( ...
    cfg.to.TOFL_ft, to.rho_slugft3, cfg.to.CLmax, to.CLg, cfg.to.CD0, to.K, ...
    cfg.atm.g0, res.eng.TW_allow_SL);

% Landing cap: Raymer landing sizing relation (simple preliminary-design form)
ld.WS_max_lbf_ft2 = WS_from_landing_Raymer_lbft2( ...
    cfg.ld.LFL_ft, ld.rho_slugft3 / cfg.atm.rho0, cfg.ld.CLmax, cfg.ld.mL_mTO);

% Cruise cap: keep cruise CL near design target
cr.q_psf = 0.5 * cr.rho_slugft3 * cr.V_fts^2;
cr.WS_max_lbf_ft2 = cr.q_psf * cfg.cr.CLdes / max(res.wt.WmidCruise_W0, 1e-6);

% Final wing-loading choice
ws.WS_ref_lbf_ft2 = min([to.WS_max_lbf_ft2, ld.WS_max_lbf_ft2, cr.WS_max_lbf_ft2]);
ws.WS_ref_Nm2 = ws.WS_ref_lbf_ft2 * 47.88025898;

%% WING AREA
ws.Sref_ft2 = res.wt.W0_lbf / ws.WS_ref_lbf_ft2;
ws.Sref_m2  = ws.Sref_ft2 * 0.09290304;

ws.b_ft     = sqrt(cfg.wg.AR * ws.Sref_ft2);
ws.cbar_ft  = ws.Sref_ft2 / ws.b_ft;

%% TAKEOFF / LANDING / CRUISE CHECKS AT FINAL WING SIZE
to.TW_req = TW_from_takeoff_Raymer( ...
    cfg.to.TOFL_ft, ws.WS_ref_lbf_ft2, to.rho_slugft3, cfg.to.CLmax, to.CLg, ...
    cfg.to.CD0, to.K, cfg.atm.g0);

to.TOFL_est_ft = takeoff_groundroll_Raymer_ft( ...
    ws.WS_ref_lbf_ft2, res.eng.TW_avail_SL, to.rho_slugft3, cfg.to.CLmax, ...
    to.CLg, cfg.to.CD0, to.K, cfg.atm.g0);

ld.WS_land_lbf_ft2 = ws.WS_ref_lbf_ft2 * cfg.ld.mL_mTO;
ld.LFL_est_ft = 80 * ld.WS_land_lbf_ft2 / max(ld.rho_slugft3 / cfg.atm.rho0 * cfg.ld.CLmax, 1e-6);

% Cruise lift / drag / thrust
cr.WS_cruise_lbf_ft2 = ws.WS_ref_lbf_ft2 * res.wt.WmidCruise_W0;
cr.CL_cruise = cr.WS_cruise_lbf_ft2 / cr.q_psf;
cr.CD_cruise = cfg.wg.CD0cr + cr.K * cr.CL_cruise^2;
cr.LD_cruise = cr.CL_cruise / cr.CD_cruise;
cr.TW_req = cr.CD_cruise / cr.CL_cruise;
cr.lapse = thrust_lapse_simple(cfg.cr.h_ft, cfg.cr.M);
cr.TW_req_SL_eq = cr.TW_req / max(cr.lapse, 1e-6);

to.pass = to.TW_req <= res.eng.TW_allow_SL;
ld.pass = ld.LFL_est_ft <= cfg.ld.LFL_ft;
cr.pass = cr.TW_req <= (res.eng.TW_allow_SL * max(cr.lapse, 1e-6));

%% SPEEDS
rho_to = to.rho_slugft3;
rho_ld = ld.rho_slugft3;

to.Vs_fts  = sqrt(2 * ws.WS_ref_lbf_ft2 / max(rho_to * cfg.to.CLmax, 1e-12));
to.V2_fts   = 1.20 * to.Vs_fts;
ld.Vs_fts   = sqrt(2 * (ws.WS_ref_lbf_ft2 * cfg.ld.mL_mTO) / max(rho_ld * cfg.ld.CLmax, 1e-12));
ld.Vapp_fts = 1.30 * ld.Vs_fts;

to.Vs_kt  = to.Vs_fts / 1.68780986;
to.V2_kt   = to.V2_fts / 1.68780986;
ld.Vs_kt   = ld.Vs_fts / 1.68780986;
ld.Vapp_kt = ld.Vapp_fts / 1.68780986;

[to.rho_chk, to.mu_chk, to.Re] = calcISA(0, ws.cbar_ft, to.V2_fts);
[ld.rho_chk, ld.mu_chk, ld.Re] = calcISA(0, ws.cbar_ft, ld.Vapp_fts);
[cr.rho_chk, cr.mu_chk, cr.Re] = calcISA(cfg.cr.h_ft, ws.cbar_ft, cr.V_fts);

%% OUTPUT STRUCT
res.to = to;
res.ld = ld;
res.cr = cr;
res.ws = ws;
res.cfg = cfg;

%% PRINT SUMMARY
fprintf('\n---- Mission / Weight ----\n');
fprintf('W0 (MTOW)                   : %.0f lbf\n', res.wt.W0_lbf);
fprintf('We (empty)                  : %.0f lbf\n', res.wt.We_lbf);
fprintf('Wf (fuel)                   : %.0f lbf\n', res.wt.Wf_lbf);
fprintf('Payload + crew              : %.0f lbf\n', res.wt.Wfix_lbf);

fprintf('\n---- Wing Size ----\n');
fprintf('Final W/S                   : %.1f lbf/ft^2\n', ws.WS_ref_lbf_ft2);
fprintf('Wing area                   : %.0f ft^2\n', ws.Sref_ft2);
fprintf('Span                        : %.1f ft\n', ws.b_ft);
fprintf('Mean chord                  : %.2f ft\n', ws.cbar_ft);

fprintf('\n---- Constraint Caps ----\n');
fprintf('Takeoff W/S cap             : %.1f lbf/ft^2\n', to.WS_max_lbf_ft2);
fprintf('Landing W/S cap             : %.1f lbf/ft^2\n', ld.WS_max_lbf_ft2);
fprintf('Cruise W/S cap              : %.1f lbf/ft^2\n', cr.WS_max_lbf_ft2);

fprintf('\n---- Thrust ----\n');
fprintf('Installed SL T/W            : %.3f\n', res.eng.TW_avail_SL);
fprintf('Allowed T/W after margin    : %.3f\n', res.eng.TW_allow_SL);
fprintf('Takeoff required T/W        : %.3f\n', to.TW_req);
fprintf('Cruise required T/W (alt)   : %.3f\n', cr.TW_req);
fprintf('Cruise req T/W SL-equiv     : %.3f\n', cr.TW_req_SL_eq);

fprintf('\n---- Field Length Checks ----\n');
fprintf('Takeoff est. ground roll    : %.0f ft\n', to.TOFL_est_ft);
fprintf('Takeoff requirement         : %.0f ft   %s\n', cfg.to.TOFL_ft, ternary(to.pass,'OK','FAIL'));
fprintf('Landing est. field length   : %.0f ft\n', ld.LFL_est_ft);
fprintf('Landing requirement         : %.0f ft   %s\n', cfg.ld.LFL_ft, ternary(ld.pass,'OK','FAIL'));
fprintf('Cruise thrust check         : %s\n', ternary(cr.pass,'OK','FAIL'));

fprintf('\n---- Speeds ----\n');
fprintf('Takeoff Vs                  : %.1f kt\n', to.Vs_kt);
fprintf('Takeoff V2                  : %.1f kt\n', to.V2_kt);
fprintf('Landing Vs                  : %.1f kt\n', ld.Vs_kt);
fprintf('Approach Vapp               : %.1f kt\n', ld.Vapp_kt);

fprintf('\n---- Reynolds Numbers ----\n');
fprintf('Takeoff Re                  : %.3e\n', to.Re);
fprintf('Landing Re                  : %.3e\n', ld.Re);
fprintf('Cruise Re                   : %.3e\n', cr.Re);


%% UI FIGURE
f = uifigure('Name','Aircraft Sizing Results');

txt = sprintf([ ...
    'MTOW: %.0f lbf\n' ...
    'Wing Area: %.0f ft^2\n' ...
    'W/S: %.1f lb/ft^2\n' ...
    'Takeoff: %.0f ft\n' ...
    'Landing: %.0f ft\n'], ...
    out.wt.W0_lbf, ...
    out.ws.Sref_ft2, ...
    out.ws.WS_ref_lbf_ft2, ...
    out.to.TOFL_est_ft, ...
    out.ld.LFL_est_ft);

uilabel(f,'Text',txt,'Position',[20 20 300 200]);




%% OPTIONAL CONSTRAINT PLOT
if cfg.opt.doPlot
    WS_grid = linspace(30, 220, 400); % lbf/ft^2

    TW_to_grid = arrayfun(@(ws) TW_from_takeoff_Raymer( ...
        cfg.to.TOFL_ft, ws, to.rho_slugft3, cfg.to.CLmax, to.CLg, ...
        cfg.to.CD0, to.K, cfg.atm.g0), WS_grid);

    TW_cr_grid = (cr.q_psf .* cfg.wg.CD0cr ./ WS_grid + cr.K .* WS_grid ./ cr.q_psf) ./ max(cr.lapse, 1e-6);
    TW_cl_grid = nan(size(WS_grid)); % add later if you want a climb line

    figure; hold on; grid on;
    plot(WS_grid, TW_to_grid, 'LineWidth', 1.8);
    plot(WS_grid, TW_cr_grid, 'LineWidth', 1.8);
    xline(ld.WS_max_lbf_ft2, '--', 'LineWidth', 1.4);
    xline(cr.WS_max_lbf_ft2, ':', 'LineWidth', 1.4);
    xline(ws.WS_ref_lbf_ft2, '-.', 'LineWidth', 1.4);
    plot(ws.WS_ref_lbf_ft2, res.eng.TW_allow_SL, 'o', 'MarkerSize', 8, 'LineWidth', 1.8);

    xlabel('W/S [lbf/ft^2]');
    ylabel('Required T/W [-]');
    title('Raymer Constraint Diagram');
    legend({'Takeoff', 'Cruise', 'Landing cap', 'Cruise cap', 'Selected W/S', 'Design point'}, ...
        'Location', 'best');
end

out = res;

end

%% ===== HELPER FUNCTIONS =====

function W0_lbf = solve_W0_lbf(Wfix_lbf, Wf_W0, A, C, W0_guess_lbf, tol)
W0_lbf = W0_guess_lbf;

for k = 1:800
    We_W0 = A * W0_lbf^C;
    denom = 1 - We_W0 - Wf_W0;

    if denom <= 0.02
        error('Infeasible weight fractions: 1 - We/W0 - Wf/W0 = %.4f', denom);
    end

    W0_new = Wfix_lbf / denom;
    err = abs(W0_new - W0_lbf) / max(W0_lbf, 1);

    W0_lbf = W0_new;

    if err < tol
        break;
    end
end
end

function WS_max_lbf_ft2 = WS_from_landing_Raymer_lbft2(LFL_ft, sigma_L, CLmax_L, mL_mTO)
% Simple Raymer landing sizing form, preliminary-design level.
% Uses landing wing loading based on landing weight, then converts to MTOW wing loading.

Sa_ft = 0;
WS_land_lbf_ft2 = max((LFL_ft - Sa_ft) * sigma_L * CLmax_L / 80.0, 0);
WS_max_lbf_ft2 = WS_land_lbf_ft2 / max(mL_mTO, 1e-6);
end

function WS_max_lbf_ft2 = WS_from_takeoff_Raymer_lbft2(TOFL_ft, rho_slugft3, CLmax_TO, CLg_TO, CD0_TO, K_TO, g_ftps2, TW_allow)
% Solve for the maximum W/S that meets TOFL with an allowed T/W.

f = @(WS) takeoff_groundroll_Raymer_ft(WS, TW_allow, rho_slugft3, CLmax_TO, CLg_TO, CD0_TO, K_TO, g_ftps2) - TOFL_ft;

WS_lo = 1.0;
WS_hi = 400.0;

f_lo = f(WS_lo);
f_hi = f(WS_hi);

nGrow = 0;
while ~(isfinite(f_lo) && isfinite(f_hi) && f_lo <= 0 && f_hi >= 0)
    if ~isfinite(f_hi) || f_hi < 0
        WS_hi = WS_hi * 1.5;
        f_hi = f(WS_hi);
    elseif f_lo > 0
        WS_lo = WS_lo / 1.5;
        f_lo = f(WS_lo);
    end

    nGrow = nGrow + 1;
    if nGrow > 80
        error('Could not bracket takeoff W/S solution.');
    end
end

for k = 1:100
    WS_mid = 0.5 * (WS_lo + WS_hi);
    f_mid = f(WS_mid);

    if ~isfinite(f_mid)
        WS_hi = WS_mid;
        continue;
    end

    if f_mid > 0
        WS_hi = WS_mid;
    else
        WS_lo = WS_mid;
    end
end

WS_max_lbf_ft2 = 0.5 * (WS_lo + WS_hi);
end

function TW_req = TW_from_takeoff_Raymer(TOFL_ft, WS_lbf_ft2, rho_slugft3, CLmax_TO, CLg_TO, CD0_TO, K_TO, g_ftps2)
% Solve for the T/W required to meet the specified TOFL at a given W/S.

f = @(TW) takeoff_groundroll_Raymer_ft(WS_lbf_ft2, TW, rho_slugft3, CLmax_TO, CLg_TO, CD0_TO, K_TO, g_ftps2) - TOFL_ft;

TW_lo = 0.02;
TW_hi = 1.50;

f_lo = f(TW_lo);
f_hi = f(TW_hi);

nGrow = 0;
while ~(isfinite(f_lo) && isfinite(f_hi) && f_lo >= 0 && f_hi <= 0)
    if ~isfinite(f_lo) || f_lo < 0
        TW_lo = max(TW_lo / 1.5, 1e-4);
        f_lo = f(TW_lo);
    elseif f_hi > 0
        TW_hi = TW_hi * 1.5;
        f_hi = f(TW_hi);
    end

    nGrow = nGrow + 1;
    if nGrow > 80
        error('Could not bracket takeoff T/W solution.');
    end
end

for k = 1:100
    TW_mid = 0.5 * (TW_lo + TW_hi);
    f_mid = f(TW_mid);

    if ~isfinite(f_mid)
        TW_lo = TW_mid;
        continue;
    end

    if f_mid > 0
        TW_lo = TW_mid;
    else
        TW_hi = TW_mid;
    end
end

TW_req = 0.5 * (TW_lo + TW_hi);
end

function Sg_ft = takeoff_groundroll_Raymer_ft(WS_lbf_ft2, TW, rho_slugft3, CLmax_TO, CLg_TO, CD0_TO, K_TO, g_ftps2)
% Raymer-style ground-roll model with the ln() form.

if WS_lbf_ft2 <= 0 || TW <= 0
    Sg_ft = Inf;
    return;
end

mu = 0.02;  % runway rolling friction, tune if you have a better value
Kt = TW - mu;

if Kt <= 0
    Sg_ft = Inf;
    return;
end

Vto_fts = 1.1 * sqrt(2 * WS_lbf_ft2 / max(rho_slugft3 * CLmax_TO, 1e-12));
Ka = (rho_slugft3 / (2 * WS_lbf_ft2)) * (mu * CLg_TO - CD0_TO - K_TO * CLg_TO^2);

if abs(Ka) < 1e-12
    Sg_ft = Vto_fts^2 / (2 * g_ftps2 * Kt);
    return;
end

arg = 1 + (Ka / Kt) * Vto_fts^2;
if arg <= 0
    Sg_ft = Inf;
    return;
end

Sg_ft = log(arg) / (2 * g_ftps2 * Ka);
end

function a_fts = isa_speed_of_sound_ftps(h_ft)
% Troposphere-only ISA speed of sound, valid through ~36,089 ft.
T0 = 518.67;      % R
L  = 0.00356616;  % R/ft
gamma = 1.4;
R = 1716.59;      % ft*lbf/(slug*R)

T = T0 - L * h_ft;
if T <= 0
    error('Altitude is outside the valid range for this ISA troposphere model.');
end

a_fts = sqrt(gamma * R * T);
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end