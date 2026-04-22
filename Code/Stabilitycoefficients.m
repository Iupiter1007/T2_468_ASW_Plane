clear; clc;

% ====================== INPUTS (FROM OPENVSP) ======================

% Wing geometry
Sw      = 1680;
b       = 140;
c_bar   = 12;
A       = 11.65;
lambda  = 0.25;

% Aerodynamic inputs from OpenVSP
CL      = input('Flight condition CL: ');
CDi     = input('Induced drag coefficient CDi: ');
cla     = input('Airfoil lift slope cla (per rad): ');
Lambda_LE = input('Leading edge sweep (deg): ')*pi/180;

Xcg     = input('CG location (% MAC 0-1): ');

% Tail inputs
Gamma_v = 45*pi/180;
Sv_panel = 160;
eta_h = 0.90;

Xach    = input('V-tail AC (% MAC): ');
Z_v     = input('Vertical tail offset (ft): ');
CLalpha_v = input('Tail lift slope (per rad): ');

Cm_alpha_fus = input('Fuselage Cm_alpha (negative): ');

% ====================== DERIVED QUANTITIES ======================

% 1. Oswald efficiency
e = CL^2 / (pi * A * CDi);

% 2. Finite wing lift slope
CLalpha_w = cla / (1 + cla/(pi*A*e));

% 3. Downwash gradient
d_epsilon_d_alpha = (2 * CLalpha_w) / (pi * A);

% 4. Quarter-chord sweep
Lambda = atan( tan(Lambda_LE) - (1/A)*((1-lambda)/(1+lambda)) );

% Wing AC
Xacw = 0.25;

% ====================== V-TAIL EFFECTIVE AREAS ======================

Sh_eff = 2 * Sv_panel * (cos(Gamma_v)^2);
Sv_eff = 2 * Sv_panel * (sin(Gamma_v)^2);

% ====================== LONGITUDINAL ======================

numerator   = CLalpha_w * Xacw - Cm_alpha_fus + ...
              eta_h * (Sh_eff/Sw) * CLalpha_v * (1 - d_epsilon_d_alpha) * Xach;

denominator = CLalpha_w + ...
              eta_h * (Sh_eff/Sw) * CLalpha_v * (1 - d_epsilon_d_alpha);

Xnp = numerator / denominator;
static_margin = Xnp - Xcg;

Cm_alpha = -denominator * (Xnp - Xcg);

% ====================== LATERAL-DIRECTIONAL ======================

Cn_beta_w = (CL^2 / (4*pi*A)) - (tan(Lambda)/(pi*A*(A+4*cos(Lambda)))) * ...
            (cos(Lambda) - A/2 - A^2/(8*cos(Lambda)) + ...
            6*(Xacw-Xcg)*sin(Lambda)/A);

Cl_beta_w_sweep = -0.0001 * A * CL;

% NOTE: replace Gamma_v with actual wing dihedral if available
Cl_beta_w_dihedral = -CLalpha_w * (5*pi/180) * (1 + 2*lambda)/(3*(1+lambda)) * 0.25;

Cl_beta_w = Cl_beta_w_sweep + Cl_beta_w_dihedral;

Cn_beta_fus = -0.03;

Cn_beta_v = CLalpha_v * (Sv_eff / Sw) * (Xach - Xcg) * sin(Gamma_v);
Cl_beta_v = -CLalpha_v * (Sv_eff / Sw) * (Z_v / b) * cos(Gamma_v);

Cn_beta = Cn_beta_w + Cn_beta_fus + Cn_beta_v;
Cl_beta = Cl_beta_w + Cl_beta_v;

% ====================== OUTPUT ======================

disp(' ');
disp('=== RESULTS (Raymer Ch. 16 - FULLY COUPLED) ===');

disp(['Oswald efficiency e    : ', num2str(e,'%6.3f')]);
disp(['CL_alpha_w             : ', num2str(CLalpha_w,'%6.3f'), ' per rad']);
disp(['Downwash dε/dα         : ', num2str(d_epsilon_d_alpha,'%6.3f')]);
disp(['Lambda_c/4 (deg)       : ', num2str(Lambda*180/pi,'%6.2f')]);

disp(['Static Margin          : ', num2str(static_margin*100, '%6.2f'), ' %']);
disp(['Cm_alpha               : ', num2str(Cm_alpha, '%6.4f'), ' per rad']);
disp(['Cn_beta                : ', num2str(Cn_beta, '%6.4f')]);
disp(['Cl_beta                : ', num2str(Cl_beta, '%6.4f')]);

% Stability checks
if static_margin > 0.05
    disp('→ Longitudinally stable');
else
    disp('→ WARNING: Low longitudinal stability');
end

if Cn_beta > 0
    disp('→ Directionally stable');
else
    disp('→ WARNING: Directional instability');
end