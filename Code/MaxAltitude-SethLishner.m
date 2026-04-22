clc; clear; close all;

%% ================== INPUT PARAMETERS ==================
W   = 180000;        % Weight (lbs)
S   = 1680;          % Wing area (ft^2)
CD0 = 0.005414*3;      % Parasite drag coefficient
k   = 0.05;          % Induced drag factor

T0   = 46883;        % Sea level thrust (lbs)
rho0 = 0.0023769;    % Sea level density (slugs/ft^3)

engine_type = 'turboprop'; % 'jet' or 'turboprop'

%% ================== ALTITUDE RANGE ==================
h_max_guess = 70000; % ft
dh = 100;            
h_vec = 0:dh:h_max_guess;

%% ================== STORAGE ==================
T_req_vec   = zeros(size(h_vec));
T_avail_vec = zeros(size(h_vec));

%% ================== LOOP ==================
for i = 1:length(h_vec)
    
    h = h_vec(i);
    
    % ---- Standard Atmosphere ----
    if h <= 36089  % Troposphere
        T = 518.67 - 0.003566*h;
        rho = rho0*(T/518.67)^((32.174/(0.003566*1716))-1);
    else           % Lower Stratosphere
        T = 389.97;
        rho = 0.000706; % approx constant density decay handled below
        
        % exponential decay above tropopause
        rho = rho * exp(-(h-36089)/20806);
    end
    
    % ---- Optimal CL (minimum thrust required) ----
    CL_opt = sqrt(CD0/k);
    
    % ---- Velocity ----
    V = sqrt( (2*W) / (rho*S*CL_opt) );
    
    % ---- Drag ----
    CD = CD0 + k*CL_opt^2;
    D  = 0.5*rho*V^2*S*CD;
    
    % ---- Thrust Required ----
    T_req = D;
    
    % ---- Thrust Available Model ----
    switch engine_type
        case 'jet'
            T_avail = T0*(rho/rho0);
        case 'turboprop'
            T_avail = T0 * (rho/rho0)^0.8 * (T/518.67)^0.5;
    end
    
    % ---- Store ----
    T_req_vec(i)   = T_req;
    T_avail_vec(i) = T_avail;
end

%% ================== FIND CEILING ==================
diff = T_avail_vec - T_req_vec;

idx = find(diff <= 0, 1, 'first');

if isempty(idx)
    disp('Increase altitude range.');
else
    % Linear interpolation for better accuracy
    h1 = h_vec(idx-1);
    h2 = h_vec(idx);
    
    d1 = diff(idx-1);
    d2 = diff(idx);
    
    h_max = h1 + (0 - d1)*(h2 - h1)/(d2 - d1);
    
    fprintf('Absolute Ceiling: %.2f ft\n', h_max);
end

%% ================== PLOT ==================
figure;
plot(h_vec, T_req_vec, 'LineWidth', 2); hold on;
plot(h_vec, T_avail_vec, '--', 'LineWidth', 2);
xlabel('Altitude (ft)');
ylabel('Thrust (lb)');
legend('Thrust Required','Thrust Available');
title('Ceiling Determination');
grid on;