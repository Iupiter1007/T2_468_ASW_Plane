% sizing_fsolve_raymer.m
% Iterative initial sizing using mission weight fractions + Raymer empty-weight estimate
% Returns W0 (lb)

clearvars; close all; clc;

%% Given payload/crew and mission
W.crew   = 800;        % lb
W.payload = 3e4;       % lb

% Fuel consumption (specific fuel consumption)
cruise.C = 0.5;    % lb fuel / (hr * lb thrust-equivalent)  -- user provided unit, used as lb/lb/hr
loiter.C = 0.4;    % lb/lb/hr for loiter (as provided)

% Aerodynamics
S_wet_S_ref = 4.5;           % wetted area ratio (unused here but kept)
L_D.max  = 19;             % maximum L/D (clean)
L_D.cruise = 0.866 * L_D.max;  % cruise L/D (given)

% Mission profile
mission.total_range_mi = 4000;   % total mission distance (mi)
range_oneway_mi = mission.total_range_mi / 2;  % one leg (mi)
cruise.range_mi = range_oneway_mi;     

% Cruise speed estimate (user expression simplified and converted to mph)
% original expression: sqrt(1.4*287*227) * 0.85 * 2.23  (units were inconsistent)
% Here we keep the user's numeric result approach but label units clearly:
cruise.vel_mph = sqrt(1.4 * 287 * 227) * 0.85 * 2.23;  % numeric mph result using same formula
% Note: if you want physically correct speed, replace with your chosen TAS in mph

loiter.time_hr = 4;    % hr

%% Non-mission fractions (typical housekeeping fractions)
% These are user-supplied or typical values:
W2_W0 = 0.97 * 0.985;   % e.g., engine start, taxi, takeoff reserves (example)
W6_W5 = 0.995;          % e.g., landing/contingency small fraction

%% Compute cruise/loiter weight fractions (ensure correct Breguet sign and consistent units)
% Range in hours for one cruise leg:
cruise.time_hr_oneway = cruise.range_mi / cruise.vel_mph;    % hr

% Cruise weight fraction for one leg (exp(-C * time / (L/D))) 
% NOTE: C must be in 1/hr (it is lb_fuel / (lb_thrust * hr) as given). This form follows the classic
% exponential fuel fraction: W_after/W_before = exp(- C * t / (L/D) )
W3_W2 = exp( - cruise.C * cruise.time_hr_oneway / L_D.cruise );  % cruise outbound
W5_W4 = W3_W2;   % return cruise assumed same

% Loiter weight fraction over loiter.time_hr (uses best L/D if available)
W4_W3 = exp( - loiter.C * loiter.time_hr / L_D.max );  % loiter stage uses L_D.max

% Compose all
W6_W0 = W2_W0 * W3_W2 * W4_W3 * W5_W4 * W6_W5;

% Fuel fraction (including some growth factor / unusable fuel factor)
Wf_W0 = 1.06 * (1 - W6_W0);   % the 1.06 factor was in user's code

%% Raymer empirical empty-weight estimate
% Raymer provides statistical fits of empty weight vs MTOW. These are empirical.
% We implement a simple power-law: We = a * W0^b  --> empty fraction = We / W0 = a * W0^(b-1)
% Default coefficients below are placeholders you can tune to match the aircraft class.
% military bomber
raymer.a = 0.93;   % empirical multiplier (tunable)
raymer.b = -0.07;   % empirical exponent (tunable)

% jet transport
% raymer.a = 1.02;   % empirical multiplier (tunable)
% raymer.b = -0.06;   % empirical exponent (tunable)

% We'll solve for W0 via fsolve on the equation:
% F(W0) = W0 - (W.crew + W.payload) / (1 - Wf_W0 - We/W0) = 0
% where We/W0 = raymer.a * W0^(raymer.b - 1)

%% Root find with fsolve
% Initial guess (lb)
W0_guess = 200000;

fun = @(W0) sizing_residual(W0, W, Wf_W0, raymer);

[W0_sol, fval, exitflag, output] = fsolve(fun, W0_guess);


% Compute breakdown at solution
We_W0 = raymer.a * W0_sol^(raymer.b);    % empty weight (lb)
Wfuel = Wf_W0 * W0_sol;
W_payload_crew = W.payload + W.crew;

fprintf('\n----- Sizing result -----\n');
fprintf('MTOW W0 = %.0f lb (%.2f klb)\n', W0_sol, W0_sol/1000);
fprintf('Empty weight We = %.0f lb (%.2f %% of MTOW)\n', We_W0*W0_sol, 100*We_W0);
fprintf('Fuel weight Wf = %.0f lb (%.2f %% of MTOW)\n', Wfuel, 100*Wf_W0);
fprintf('Payload+Crew = %.0f lb (%.2f %% of MTOW)\n', W_payload_crew, 100*W_payload_crew/W0_sol);
fprintf('Check: fractions sum = %.4f (should be ~1)\n', We_W0 + Wf_W0 + W_payload_crew / W0_sol);
fprintf('Mission fractions: W2/W0=%.4f, cruise leg one=%.4f, loiter=%.4f\n', W2_W0, W3_W2, W4_W3);

%% Local function
function F = sizing_residual(W0, W, Wf_W0, raymer)
    % W0 scalar (lb)
    % W struct with fields crew, payload
    % Wf_W0 scalar (fuel fraction of W0 determined by mission)
    % raymer struct with a,b such that We = a * W0^b
    
    if W0 <= 0
        F = -W0; return;
    end
    
    We_W0 = raymer.a * W0^(raymer.b);      % empty weight (lb)
    
    % Equation to solve: W0 = (W.crew + W.payload) / (1 - Wf_W0 - We_W0)
    % Rearranged residual:
    F = W0 - (W.crew + W.payload) / (1 - Wf_W0 - We_W0);
end
