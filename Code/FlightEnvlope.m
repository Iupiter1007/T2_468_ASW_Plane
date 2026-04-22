clc; clear; close all;

%% ===================== INPUTS =====================

% Payload (lb)
W_ISR     = 10000;
W_mun     = 20000;
W_crew    = 800;
W_payload = W_ISR + W_mun;

% Mission
Range_leg_nmi = 2000;
E_loiter_hr   = 4.0;
E_reserve_hr  = 0.75;
fuel_cont     = 0.05;

% Cruise
Mach_cruise = 0.75;
h_cruise_ft = 30000;

% Aero
LD_cruise = 18;
LD_loiter = 16;
c_cruise  = 0.55;
c_loiter  = 0.50;

% Field requirements
TOFL_req_ft = 6000;
LFL_req_ft  = 5000;

CLmax_TO = 2.5;
CLmax_L  = 2.8;
mML_over_mMTO = 0.78;

sigma_TO = 1.0;
sigma_L  = 1.0;

% Empty weight regression
A = 0.93; C = -0.07;

% Engine
LEAP_rating_lbf_SL = 31000;
nEng = 2;

thrust_margin = 1.20;
margin_WS     = 0.95;

CL_design_SC2 = 0.40;

W0_guess = 170000;
tol      = 1e-4;

rhoSL     = 1.225;
kt_per_ms = 1.94384;
g0        = 9.80665;

kL  = 0.107;
kTO = 2.34;

%% ===================== OUTER LOOP =====================

for outer = 1:40

    V_cruise_kt = mach_to_tas_knots(Mach_cruise, h_cruise_ft);
    t_leg_hr    = Range_leg_nmi / V_cruise_kt;

    W1_W0=0.970; W2_W1=0.985; W3_W2=0.995;
    W7_W6=0.995; W8_W7=0.980;

    W4_W3 = exp(-(t_leg_hr * c_cruise) / LD_cruise);
    W5_W4 = exp(-(E_loiter_hr * c_loiter) / LD_loiter);
    W6_W5 = exp(-(t_leg_hr * c_cruise) / LD_cruise);
    W9_W8 = exp(-(E_reserve_hr * c_loiter) / LD_loiter);

    Mff_nom = W1_W0*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5*W7_W6*W8_W7*W9_W8;
    Wf_W0_nom = 1 - Mff_nom;
    Wf_W0 = 1 - (1 - Wf_W0_nom)*(1 - fuel_cont);

    W3_W0 = W1_W0*W2_W1*W3_W2;
    WmidCruise_over_W0 = W3_W0*sqrt(W4_W3);

    sigma_cruise = density_ratio_ISA(h_cruise_ft);
    rho_cruise   = rhoSL * sigma_cruise;

    V_cruise_ms = V_cruise_kt/kt_per_ms;
    q_cruise = 0.5*rho_cruise*V_cruise_ms^2;

    WS_L = WS_from_landing(LFL_req_ft,sigma_L,CLmax_L,mML_over_mMTO,kL);
    WS_C = (q_cruise*CL_design_SC2)/WmidCruise_over_W0;

    WS_MTO = margin_WS*min(WS_L,WS_C);

    mS_MTO = WS_MTO/g0;

    TW_req = mS_MTO*(kTO*sigma_TO)/(ft_to_m(TOFL_req_ft)*CLmax_TO);
    TW_design = thrust_margin*TW_req;

    W0 = solve_W0(W_payload+W_crew,Wf_W0,A,C,W0_guess,tol);

    W0_N = lb_to_N(W0);
    S_m2 = W0_N/WS_MTO;
    S_ft2 = m2_to_ft2(S_m2);

    T_total = TW_design*W0;
    T_each  = T_total/nEng;

    if abs(W0-W0_guess)/W0_guess < 5e-4
        break
    end
    W0_guess=W0;
end

%% ===================== EXPANDED FLIGHT ENVELOPE =====================
figure; hold on; grid on; box on;

n_limit_pos = 2.5;
n_limit_neg = -1.0;

h_env_ft = 0;
sigma_env = density_ratio_ISA(h_env_ft);
rho_env = rhoSL*sigma_env;
rho_slug = rho_env*0.00194032;

CLmax_clean = 1.6;
CL_alpha = 5.5;
Ude = 56;

V_C_kt = mach_to_tas_knots(Mach_cruise,h_env_ft);
V_D_kt = 1.25*V_C_kt;
V_kt = linspace(0,V_D_kt,800);
V_ft_s = V_kt/0.592484;

%% STALL
n_stall = (rho_slug.*V_ft_s.^2*S_ft2*CLmax_clean)/(2*W0);
n_pos = min(n_stall,n_limit_pos);
n_neg = -min(n_stall,abs(n_limit_neg));

%% GUST
mu = (2*W0)/(rho_slug*S_ft2*CL_alpha*32.2);
Kg = (0.88*mu)/(5.3+mu);

n_gust = 1 + (Kg*rho_slug*V_ft_s*Ude*CL_alpha*S_ft2)/(2*W0);
n_gust_neg = 1 - (Kg*rho_slug*V_ft_s*Ude*CL_alpha*S_ft2)/(2*W0);

%% OPERATIONAL AREA
n_upper = min(n_pos, n_gust);
n_lower = max(n_neg, n_gust_neg);

mask_valid = n_upper > n_lower;

fill([V_kt(mask_valid) fliplr(V_kt(mask_valid))], ...
     [n_upper(mask_valid) fliplr(n_lower(mask_valid))], ...
     [0 0.7 0], 'FaceAlpha',0.3, 'EdgeColor','none');

%% THRUST LIMIT
T_available = T_total;

CD0 = 0.02;
AR = 9;
e = 0.8;

q = 0.5*rho_slug.*V_ft_s.^2;
CL = W0./(q*S_ft2);
CD = CD0 + (CL.^2)/(pi*AR*e);

D = q*S_ft2.*CD;

V_thrust_limit = V_kt(D > T_available);

if ~isempty(V_thrust_limit)
    V_min_thrust = min(V_thrust_limit);
else
    V_min_thrust = 0;
end

%% ===================== FIXED DYNAMIC PRESSURE LIMIT =====================
V_D_ft_s = V_D_kt / 0.592484;
q_limit = 0.5 * rho_slug * V_D_ft_s^2;
V_q_limit = V_D_kt;

%% PLOTS
plot(V_kt,n_pos,'b','LineWidth',2);
plot(V_kt,n_neg,'b','LineWidth',2);

plot(V_kt,n_gust,'m--','LineWidth',1.5);
plot(V_kt,n_gust_neg,'m--','LineWidth',1.5);

yline(n_limit_pos,'r--','LineWidth',1.5);
yline(n_limit_neg,'r--','LineWidth',1.5);

xline(V_C_kt,'k-.','LineWidth',1.5);
xline(V_D_kt,'k','LineWidth',2);

if V_min_thrust > 0
    xline(V_min_thrust,'g--','LineWidth',1.5);
end

xline(V_q_limit,'c--','LineWidth',1.5);

xlabel('Velocity (knots)');
ylabel('Load Factor (n)');
title('Expanded Flight Envelope (Including Gust, Thrust & q Limits)');

legend('Operational Area (Gust + Stall)',...
       'Positive Stall','Negative Stall',...
       'Gust +','Gust -',...
       'n_{limit}','n_{neg}',...
       'V_C','V_D',...
       'Thrust Limit','q Limit',...
       'Location','bestoutside');

ylim([n_limit_neg-1 n_limit_pos+1]);

fprintf('\n============= EXPANDED FLIGHT ENVELOPE =============\n');
fprintf('Cruise speed Vc: %.1f kt\n',V_C_kt);
fprintf('Dive speed Vd: %.1f kt\n',V_D_kt);
fprintf('Dynamic pressure limit speed: %.1f kt\n',V_q_limit);
fprintf('Thrust-limited minimum speed: %.1f kt\n',V_min_thrust);

%% FUNCTIONS
function W0=solve_W0(Wfixed,Wf_W0,A,C,W0_init,tol)
W0=W0_init; err=1;
while err>tol
We_W0=A*W0^C;
W0_new=Wfixed/(1-We_W0-Wf_W0);
err=abs(W0_new-W0)/W0;
W0=W0_new;
end
end

function WS=WS_from_landing(LFL_ft,sigma,CLmax,mfrac,kL)
g0=9.80665;
mS=kL*sigma*CLmax*ft_to_m(LFL_ft);
WS=(mS/mfrac)*g0;
end

function sigma=density_ratio_ISA(h_ft)
rho0=1.225; T0=288.15; p0=101325;
g0=9.80665; R=287.05; a=-0.0065;
h=h_ft*0.3048;
if h<=11000
T=T0+a*h; p=p0*(T/T0)^(-g0/(a*R));
else
T11=T0+a*11000;
p11=p0*(T11/T0)^(-g0/(a*R));
T=T11; p=p11*exp(-g0*(h-11000)/(R*T));
end
rho=p/(R*T);
sigma=rho/rho0;
end

function V=mach_to_tas_knots(M,h_ft)
gamma=1.4; R=287.05; T0=288.15; a=-0.0065;
h=h_ft*0.3048;
if h<=11000
T=T0+a*h;
else
T=T0+a*11000;
end
a_sound=sqrt(gamma*R*T);
V=M*a_sound*1.94384;
end

function m=ft_to_m(ft), m=ft*0.3048; end
function ft2=m2_to_ft2(m2), ft2=m2/(0.3048^2); end
function N=lb_to_N(lb), N=lb*4.44822; end