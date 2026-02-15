% Team Leviathan
% Iterative initial sizing using mission weight fractions + Raymer empty-weight estimate
% Returns W0 (lb)

clearvars; close all; clc;

%% Given specs
W.crew   = 800;        % lb
W.payload = 3e4;       % lb

% Mission profile
mission.total_range_mi = 4000;   % total mission distance (mi)
cruise.range_mi = mission.total_range_mi / 2;  % one leg (mi)
loiter.time_hr = 4;    % hr

%% Non-mission fractions (typical housekeeping fractions)

% Fuel consumption (specific fuel consumption)
cruise.C = 0.5;    % lb fuel / (hr * lb thrust-equivalent)  
loiter.C = 0.4;    % lb/lb/hr for loiter

% Aerodynamics
S_wet_S_ref = 4.5;           % wetted area ratio (unused here but kept)
L_D.max  = 19;             % maximum L/D (clean)
L_D.cruise = 0.866 * L_D.max;  % cruise L/D (given)
% Cruise speed estimate assume temperature at 30k ft
cruise.vel_mph = sqrt(1.4 * 287 * 227) * 0.85 * 2.23;


%% Weight fractions
W2_W0 = 0.97 * 0.985;   % e.g., engine start, taxi, takeoff reserves (example)
W6_W5 = 0.995;          % e.g., landing/contingency small fraction

% Cruise weight fraction for one leg
W3_W2 = exp( - cruise.C * cruise.range_mi / cruise.vel_mph / L_D.cruise );  % cruise outbound
W5_W4 = W3_W2;   % return cruise assumed same

W4_W3 = exp( - loiter.C * loiter.time_hr / L_D.max );  % loiter stage uses L_D.max

W6_W0 = W2_W0 * W3_W2 * W4_W3 * W5_W4 * W6_W5;

Wf_W0 = 1.06 * (1 - W6_W0);  % Fuel fraction

%% Raymer empirical empty-weight estimate
% Raymer provides statistical fits of empty weight vs MTOW. These are empirical.
% We/W0 = a * W0^b
% Default coefficients below are placeholders you can tune to match the aircraft class.
% military bomber
raymer.a = 0.97;   % empirical multiplier (tunable)
raymer.b = -0.07;   % empirical exponent (tunable)

% jet transport
% raymer.a = 1.02;   % empirical multiplier (tunable)
% raymer.b = -0.06;   % empirical exponent (tunable)

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
    We_W0 = raymer.a * W0^(raymer.b);      % empty weight (lb)
    F = W0 - (W.crew + W.payload) / (1 - Wf_W0 - We_W0);
end
