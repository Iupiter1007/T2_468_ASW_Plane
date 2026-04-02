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
