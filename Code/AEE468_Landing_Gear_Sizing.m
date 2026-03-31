% AEE468_Landing_Gear_Sizing.m
% Raymer-style conceptual landing gear sizing
clear; clc;

g = 9.80665; % m/s^2

fprintf('=== Raymer-style Conceptual Landing Gear Sizing ===\n');
fprintf('All lengths in meters (m), masses in kg\n\n');

%% ------------------ USER INPUTS ------------------

MTOW_kg = input_with_default('Enter MTOW (kg) [18000]: ',18000);
W_landing_kg = input_with_default('Enter landing mass (kg) [0.94*MTOW]: ',0.94*MTOW_kg);

x_CG = input_with_default('Enter aft CG from nose (m) [8.5]: ',8.5);
x_NG = input_with_default('Enter nose gear location from nose (m) [1.0]: ',1.0);

x_tail = input_with_default('Enter tail-most point location from nose (m) [14.5]: ',14.5);
h_tail = input_with_default('Enter tail height above ground (m) [0.9]: ',0.9);
h_CG = input_with_default('Enter CG height above ground (m) [1.5]: ',1.5);

f_nose = input_with_default('Desired nose load fraction (0.08-0.15) [0.10]: ',0.10);
phi_overturn_deg = input_with_default('Target overturn angle (deg) [30]: ',30);
V_sink = input_with_default('Sink rate at touchdown (m/s) [3]: ',3);
eta = input_with_default('Shock efficiency (0.7-0.9) [0.85]: ',0.85);

n_mains = input_with_default('Number of main gear legs [2]: ',2);
wheels_per_main = input_with_default('Wheels per main leg [2]: ',2);
LLF = input_with_default('Landing load factor (1.5-2.0) [1.5]: ',1.5);

gear_weight_fraction = input_with_default('Gear weight fraction of MTOW (~0.04) [0.04]: ',0.04);
theta_rot_avail_deg = input_with_default('Available rotation angle (deg) [12]: ',12);

%% ------------------ CALCULATIONS ------------------

W_total = W_landing_kg * g;
f_mains = 1 - f_nose;

% Solve for main gear longitudinal location
x_MG = (x_CG - f_nose * x_NG) / (1 - f_nose);
wheelbase = x_NG - x_MG;

% Load split
W_mains_total = f_mains * W_total;
W_nose_total = f_nose * W_total;

n_main_wheels = n_mains * wheels_per_main;

W_per_main_static = W_mains_total / n_main_wheels;
W_per_main_peak = W_per_main_static * LLF;

% Overturn track requirement
phi_rad = deg2rad(phi_overturn_deg);
track_required = 2 * h_CG * tan(phi_rad);

% Tipback / tailstrike
theta_tipback_deg = rad2deg( atan( h_tail / (x_tail - x_MG) ) );
tailstrike = theta_tipback_deg <= theta_rot_avail_deg;

% Shock stroke (energy method)
s_stroke = 0.5 * V_sink^2 / (f_mains * g * eta);
s_in = s_stroke * 39.37;

% Gear weight estimate
gear_weight_kg = gear_weight_fraction * MTOW_kg;

%% ------------------ OUTPUT ------------------

fprintf('\n===== RESULTS =====\n');
fprintf('Main Gear Location (x_MG): %.3f m\n', x_MG);
fprintf('Wheelbase: %.3f m\n', wheelbase);

fprintf('\nLoad Distribution:\n');
fprintf('  Nose Load: %.1f %%\n', f_nose*100);
fprintf('  Per Main Wheel Static Load: %.1f N\n', W_per_main_static);
fprintf('  Per Main Wheel Peak Load: %.1f N\n', W_per_main_peak);

fprintf('\nLateral Stability:\n');
fprintf('  Required Track Width: %.3f m\n', track_required);

fprintf('\nTipback Angle: %.2f deg\n', theta_tipback_deg);
if tailstrike
    fprintf('  WARNING: Tailstrike risk at chosen rotation angle.\n');
else
    fprintf('  Rotation angle acceptable.\n');
end

fprintf('\nShock Absorber Stroke Required:\n');
fprintf('  %.3f m (%.1f inches)\n', s_stroke, s_in);

fprintf('\nEstimated Gear Weight:\n');
fprintf('  %.1f kg\n', gear_weight_kg);

fprintf('\n--- Use tire manufacturer charts with peak wheel loads for final tire sizing ---\n');

%% ------------------ HELPER FUNCTION ------------------

function val = input_with_default(prompt, default)
    user_input = input(prompt,'s');
    if isempty(user_input)
        val = default;
    else
        val = str2double(user_input);
        if isnan(val)
            val = default;
        end
    end
end