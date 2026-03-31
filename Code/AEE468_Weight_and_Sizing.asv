% MTOW_estimation_fixed_ASW.m
% Revised conceptual MTOW estimation for an ASW/ISR aircraft (Raymer-style)
% Fixes and improvements to the original script:
%  - Clarifies one-way vs round-trip cruise distance (R_oneway)
%  - Computes mission segments sequentially (outbound cruise, return cruise, loiter,
%    climb/taxi/reserves) rather than lumping ambiguous distances
%  - Adds baseline (less-conservative) scenario for quick sanity checks
%  - Defensive checks and informative diagnostics when denominator <= 0
%  - Prints intermediate fuel fractions and selected reference scenarios
%  - Keeps the L/D and empty-weight fraction trade plots
%
% Author: Steve Von Brahn (revised)
% Date: 2026-02-16 (revised)
clear; close all; clc;

%% ======= USER INPUTS =======================================
% Payload components (lb)
W_ISR   = 10000;      % ISR electronics weight [lb]
W_weap  = 20000;      % Munitions weight [lb]
W_crew  = 800;        % Crew/avionics equivalent [lb] (use 0 if unmanned)
W_payload = W_ISR + W_weap + W_crew;  % total payload [lb]

% Mission profile (clear one-way vs round-trip)
R_oneway = 2000;      % one-way cruise distance [nmi]
doRoundTrip = true;   % if true, script assumes outbound+return cruise

% Loiter and speeds
E_loiter = 4.0;       % loiter time [hr]
V_cruise = 420;       % cruise speed [ktas]

% Engine / aerodynamic guesses
c_cruise = 0.50;      % TSFC in 1/hr (Breguet units)

% Non-cruise segment fractions (product of climb/taxi/landing/reserves)
% Instead of a single guess use a small segment breakdown and product them
seg_taxi = 0.997;     % taxi fraction (W_after / W_before)
seg_climb = 0.990;    % climb + acceleration
seg_loiter_reserve = 0.995; % approach/holding reserve multiplier
non_cruise_fraction_default = seg_taxi * seg_climb * seg_loiter_reserve; % ~0.982

% Option: allow user to override with a single multiplier
use_single_non_cruise = false;  % set true to use non_cruise_fraction_single
non_cruise_fraction_single = 0.99; % if used, this will replace the product above

if use_single_non_cruise
    non_cruise_fraction = non_cruise_fraction_single;
else
    non_cruise_fraction = non_cruise_fraction_default;
end

% Sweep ranges (design options)
LoverD_vec = 14:0.25:18;              % L/D sweep
fE_vec = 0.38:0.002:0.56;              % empty-weight fraction sweep (includes 0.402 precisely)

% Fixed plotting choices
showContourLevels = 12;

%% ======= Derived mission values ==================================
if doRoundTrip
    nCruiseLegs = 2; % outbound + return
else
    nCruiseLegs = 1; % only one cruise leg
end

%% Preallocate results matrix
MTOW = nan(length(fE_vec), length(LoverD_vec));  % rows: fE, cols: L/D
FuelFraction = nan(size(MTOW));

for i = 1:length(fE_vec)
    fE = fE_vec(i);  % empty weight fraction guess
    for j = 1:length(LoverD_vec)
        LbyD = LoverD_vec(j);

        % Sequential mission: multiply weight ratios per segment
        Wratio = 1.0;

        % Cruise legs sequentially (handles 1 or 2 legs clearly)
        for leg = 1:nCruiseLegs
            cruise_exp_leg = (R_oneway * c_cruise) / (V_cruise * LbyD);
            Wratio = Wratio * exp(-cruise_exp_leg);
        end

        % Loiter (endurance form)
        loiter_exp = (c_cruise * E_loiter) / LbyD;
        Wratio = Wratio * exp(-loiter_exp);

        % Non-cruise multiplier (taxi/climb/descent/reserves)
        Wratio = Wratio * non_cruise_fraction;

        % Fuel fraction consumed during mission
        fFuel = 1 - Wratio;

        % Defensive check and MTOW solution
        denom = 1 - fE - fFuel;
        if denom <= 0
            MTOW(i,j) = NaN;   % infeasible for these assumptions
            FuelFraction(i,j) = fFuel;
        else
            MTOW(i,j) = W_payload / denom;
            FuelFraction(i,j) = fFuel;
        end
    end
end

%% ======= Plotting: Contour of MTOW over (L/D, fE)
[LbyD_grid, fE_grid] = meshgrid(LoverD_vec, fE_vec);

figure('Units','normalized','Position',[0.08 0.1 0.8 0.7]);
contourf(LbyD_grid, fE_grid, MTOW/1000, showContourLevels,'LineColor','none');  % MTOW in 1000 lb
cb = colorbar;
cb.Label.String = 'MTOW (1000 lb)';
cb.Label.FontSize = 12;

xlabel('Cruise L/D (assumed)','FontSize',12);
ylabel('Empty-weight fraction f_E = W_E / W_0','FontSize',12);
title('Estimated MTOW Trade Space','FontSize',13);
set(gca,'FontSize',11);

% Overlay a few sample points and annotate a "reference" design (choose sensible defaults)
hold on;
ref_LbyD = 18; ref_fE = 0.402; % recommended reference
ref_MTOW = interp2(LbyD_grid, fE_grid, MTOW, ref_LbyD, ref_fE);
if ~isnan(ref_MTOW)
    % Plot reference marker
plot(ref_LbyD, ref_fE,'ko', ...
    'MarkerFaceColor','w', ...
    'MarkerSize',8, ...
    'LineWidth',1.5);

% Create label string
labelStr = sprintf('Reference Design\nL/D = %.1f\nf_E = %.3f\nMTOW = %.0f lb', ...
    ref_LbyD, ref_fE, ref_MTOW);

% Place text to the LEFT of the point
text(ref_LbyD - 0.3, ref_fE, labelStr, ...
    'HorizontalAlignment','right', ...   % <-- forces text to extend left
    'VerticalAlignment','middle', ...
    'FontSize',10, ...
    'FontWeight','bold', ...
    'BackgroundColor','w', ...
    'EdgeColor','k');
end
hold off;

%% ======= Plotting: L/D sensitivity lines for a few fE values ============
figure('Units','normalized','Position',[0.12 0.12 0.75 0.6]);
hold on;
plotStyles = {'-','--','-.',':'};
fE_samples = [0.402, 0.48, 0.50, 0.54];
for k = 1:length(fE_samples)
    [~, idx] = min(abs(fE_vec - fE_samples(k)));
    plot(LoverD_vec, MTOW(idx,:)/1000, 'LineWidth',1.8,'DisplayName',sprintf('f_E=%.2f',fE_samples(k)));
end
xlabel('Cruise L/D');
ylabel('MTOW (1000 lb)');
title('MTOW vs L/D for selected empty-weight fractions');
legend('Location','northeast');
grid on; set(gca,'FontSize',11); hold off;

%% ======= Display a small results table for a chosen param set =============
fprintf('--- Example reference results (payload = %d lb) ---\n', W_payload);
fprintf('Assumptions: V=%.0f kt, R_oneway=%.0f nmi (round-trip: %s), Loiter=%.1f hr, TSFC=%.3f 1/hr\n', ...
    V_cruise, R_oneway, string(doRoundTrip), E_loiter, c_cruise);
fprintf('Non-cruise fraction used = %.4f (product of segments or single override)\n', non_cruise_fraction);
fprintf('Reference L/D=%.1f, fE=%.2f -> MTOW ~ %.0f lb (if feasible)\n', ref_LbyD, ref_fE, ref_MTOW);

% Show a small matrix for L/D = 14, 16, 18 and fE = 0.45, 0.50, 0.54
report_LDs = [14, 16, 18];
report_fEs = [0.402, 0.50, 0.54];
disp('MTOW table (lb): rows = fE, cols = L/D');
mtab = nan(length(report_fEs), length(report_LDs));
for a = 1:length(report_fEs)
    for b = 1:length(report_LDs)
        i = find(abs(fE_vec - report_fEs(a))<1e-6,1);
        j = find(abs(LoverD_vec - report_LDs(b))<1e-6,1);
        if isempty(i) || isempty(j)
            mtab(a,b) = NaN;
        else
            mtab(a,b) = MTOW(i,j);
        end
    end
end

disp(array2table(mtab,'VariableNames',strcat('L_D_',string(report_LDs)),'RowNames',strcat('fE_',string(report_fEs))));

%% ======= Print diagnostic warnings for infeasible regions ==================
% Find any infeasible combinations and display a summary
[fi,fj] = find(isnan(MTOW));
if ~isempty(fi)
    unique_fE_bad = unique(fE_vec(fi));
    fprintf('\nWARNING: infeasible region(s) detected where 1 - fE - fFuel <= 0.\n');
    fprintf('Example problematic empty-weight fractions (fE): ');
    fprintf('%.3f ', unique_fE_bad); fprintf('\n');
    fprintf('Consider reducing payload, reducing fE, increasing L/D, or reducing required range/loiter.\n');
end

%% ======= Baseline (less-conservative) scenario for quick sanity check ======
% A quick compare block using more optimistic assumptions
baseline.fE = 0.45;
baseline.R_oneway = 2000;
baseline.doRoundTrip = true;
baseline.non_cruise_fraction = 0.99;
baseline.c_cruise = c_cruise; baseline.V_cruise = V_cruise; baseline.E_loiter = E_loiter;

% Compute baseline fuel fraction and MTOW for reference L/D
LbyD_ref = ref_LbyD;
Wratio_baseline = 1.0;
for leg = 1:(baseline.doRoundTrip*1 + (~baseline.doRoundTrip)*0)
    Wratio_baseline = Wratio_baseline * exp(- (baseline.R_oneway * baseline.c_cruise) / (baseline.V_cruise * LbyD_ref));
end
Wratio_baseline = Wratio_baseline * exp(- (baseline.c_cruise * baseline.E_loiter) / LbyD_ref);
Wratio_baseline = Wratio_baseline * baseline.non_cruise_fraction;
fFuel_baseline = 1 - Wratio_baseline;
if 1 - baseline.fE - fFuel_baseline > 0
    MTOW_baseline = W_payload / (1 - baseline.fE - fFuel_baseline);
else
    MTOW_baseline = NaN;
end
fprintf('\nBaseline check (fE=%.2f, non_cruise=%.2f): fFuel=%.3f, MTOW=%.0f lb\n', baseline.fE, baseline.non_cruise_fraction, fFuel_baseline, MTOW_baseline);

%% ======= Notes & Next steps printed to workspace =======================
fprintf('\nNOTES:\n');
fprintf('- Confirm whether R_oneway is one-way or if mission requires additional legs (SAR, loiter during transit, etc.).\n');
fprintf('- Validate fE with structural/weight-statistic methods or choose reference aircraft to derive fE.\n');
fprintf('- TSFC is sensitive. Use cruise and loiter SFC appropriate for your chosen engine and altitude.\n');
fprintf('- For higher fidelity, break climb/takeoff/landing into their Raymer segment fractions (separate climb/accel, cruise-out, descent, loiter, reserve).\n');

% End of script
