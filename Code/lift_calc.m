clc; clear; close all;
%% ===================== INPUT =====================
S    = 1680;      % wing area [ft^2]
span = 140;       % span [ft]
rt_c = 12;        % reference length [ft]
W    = 180000;    % weight [lbf]

% AoA vs Cl from OpenVSP (cruise config, no flaps)
AoA = [-2.000000000000
-1.000000000000
 0.000000000000
 1.000000000000
 2.000000000000
 3.000000000000
 4.000000000000
 5.000000000000
 6.000000000000]';
Cl  = [00.051939826409
0.156319558669
0.258904087430
0.363989186245
0.475688526768
0.589465991180
0.705578572262
0.816814514041
0.855156210558]';   % <-- YOU EDIT THIS
% Speeds
cr.v_mph = 0.7*1125;
to.v_mph = 166;

%% ===================== CRUISE =====================
[cr.rho, cr.mu, cr.Re] = calcISA(30000, rt_c, cr.v_mph);
cr.v_fts = cr.v_mph * 1.4666666667;
cr.q = 0.5 * cr.rho * cr.v_fts^2;

% Lift for each AoA
cr.L   = cr.q * S .* Cl;
cr.L_W = cr.L ./ W;

% Find AoA where L/W = 1
cr.AoA_req = interp1(cr.L_W, AoA, 1);

%% ===================== TAKEOFF =====================
[to.rho, to.mu, to.Re] = calcISA(0, rt_c, to.v_mph);
to.v_fts = to.v_mph * 1.4666666667;
to.q = 0.5 * to.rho * to.v_fts^2;

% Use cruise AoA
to.Cl = interp1(AoA, Cl, cr.AoA_req);

to.L   = to.q * S * to.Cl;
to.L_W = to.L / W;

% Required Cl for takeoff (L = W)
to.Cl_req = W / (to.q * S);
to.AoA_req = interp1(Cl, AoA, to.Cl_req);

%% ===================== PRINT =====================
fprintf('\n=============== AIRCRAFT REFERENCE ===============\n');
fprintf('Wing Area (S)        : %.1f ft^2\n', S);
fprintf('Span (b)             : %.1f ft\n', span);
fprintf('Reference Length     : %.2f ft\n', rt_c);
fprintf('Weight (W)           : %.0f lbf\n', W);

fprintf('\n=============== CRUISE CONDITION ================\n');
fprintf('Velocity             : %.1f mph\n', cr.v_mph);
fprintf('Density              : %.5f slug/ft^3\n', cr.rho);
fprintf('Reynolds Number      : %.3e\n', cr.Re);
fprintf('Required AoA (L/W=1) : %.2f deg\n', cr.AoA_req);

% Table of AoA sweep
disp(table(AoA', Cl', cr.L', cr.L_W', ...
    'VariableNames', {'AoA_deg','Cl','Lift_lbf','L_W'}));

fprintf('\n=============== TAKEOFF CONDITION ===============\n');
fprintf('Velocity             : %.1f mph\n', to.v_mph);
fprintf('Density              : %.5f slug/ft^3\n', to.rho);
fprintf('Reynolds Number      : %.3e\n', to.Re);
fprintf('Using Cruise AoA     : %.2f deg\n', cr.AoA_req);
fprintf('Cl at that AoA       : %.3f\n', to.Cl);
fprintf('Lift                 : %.0f lbf\n', to.L);
fprintf('L/W                  : %.3f\n', to.L_W);

fprintf('\n--- Required for Takeoff (L=W) ---\n');
fprintf('Required Cl          : %.3f\n', to.Cl_req);
fprintf('Required AoA         : %.2f deg\n', to.AoA_req);

fprintf('\n=================================================\n');

function [rho, mu, Re] = calcISA(h_ft, Lref_ft, V_fts)
% calcISA  ISA troposphere in imperial units
%
% Inputs
%   h_ft    - altitude [ft]
%   Lref_ft - reference length [ft] (optional, for Reynolds number)
%   V_fts   - speed [fts] (optional, for Reynolds number)
%
% Outputs
%   rho - air density [slug/ft^3]
%   mu  - dynamic viscosity [slug/(ft*s)]
%   Re  - Reynolds number [-] (NaN if not requested)
%
% Notes
%   This is the standard troposphere model, good to about 36,089 ft.

    % ISA constants (imperial)
    T0 = 518.67;      % sea level temperature [R]
    L  = 0.00356616;  % lapse rate [R/ft]
    P0 = 2116.22;     % sea level pressure [lbf/ft^2]
    g  = 32.174;      % gravity [ft/s^2]
    R  = 1716.59;     % gas constant for dry air [ft*lbf/(slug*R)]

    % Dynamic viscosity constants (Sutherland's law)
    mu0 = 3.737e-7;   % reference viscosity [slug/(ft*s)] at T0
    S   = 198.72;     % Sutherland constant [R]

    % Temperature at altitude
    T = T0 - L*h_ft;

    if T <= 0
        error('Altitude is outside the valid range for this simple ISA model.');
    end

    % Pressure and density
    P = P0 * (T/T0)^(g/(L*R));
    rho = P / (R*T);

    % Dynamic viscosity
    mu = mu0 * (T/T0)^(3/2) * (T0 + S) / (T + S);

    % Optional Reynolds number
    Re = NaN;
    if nargin >= 3
        Re = rho * V_fts * Lref_ft / mu;
    end
end
