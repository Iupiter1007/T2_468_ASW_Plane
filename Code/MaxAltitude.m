clc; clear; close all;

%% ================== INPUT PARAMETERS ==================
W   = 180000;        % Weight (lbs)
S   = 1680;          % Wing area (ft^2)
CD0 = 0.005414;      % Parasite drag coefficient
k   = 0.05;          % Induced drag factor

T0  = 46883;          % Thrust available at sea level (lbs)
rho0 = 0.0023769;     % Air density at sea level (slugs/ft^3)

%% ================== ALTITUDE RANGE ==================
h_max_guess = 50000; % ft
dh = 50;              % step size

h_vec = 0:dh:h_max_guess;

%% ================== STORAGE ==================
T_req_vec   = zeros(size(h_vec));
T_avail_vec = zeros(size(h_vec));

%% ================== LOOP THROUGH ALTITUDE ==================
for i = 1:length(h_vec)
    
    h = h_vec(i);
    
    % ---- Standard Atmosphere (Troposphere) ----
    T = 518.67 - 0.003566*h; % Rankine
    rho = rho0*(T/518.67)^((32.174/(0.003566*1716))-1); % slugs/ft^3
    
    % ---- Optimal CL for minimum thrust ----
    CL_opt = sqrt(CD0/k);
    
    % ---- Velocity at lift equilibrium ----
    V = sqrt( (2*W) / (rho*S*CL_opt) ); % ft/s
    
    % ---- Drag Calculation ----
    CD = CD0 + k*CL_opt^2;
    D  = 0.5*rho*V^2*S*CD; % lbs
    
    % ---- Thrust Required ----
    T_req = D;
    
    % ---- Thrust Available (jet scaling) ----
    T_avail = T0*(rho/rho0);
    
    % ---- Store ----
    T_req_vec(i)   = T_req;
    T_avail_vec(i) = T_avail;
end

%% ================== FIND MAX ALTITUDE ==================
diff = T_avail_vec - T_req_vec;

idx = find(diff <= 0, 1, 'first');

if isempty(idx)
    disp('Max altitude not found in range. Try increasing h_max_guess.');
else
    h_max = h_vec(idx);
    fprintf('Maximum Altitude (Absolute Ceiling): %.2f ft\n', h_max);
end

%% ================== PLOT ==================
figure;
plot(h_vec, T_req_vec, 'LineWidth', 2); hold on;
plot(h_vec, T_avail_vec, '--', 'LineWidth', 2);
xlabel('Altitude (ft)');
ylabel('Thrust (lb)');
legend('Thrust Required','Thrust Available');
title('Max Altitude Determination (Imperial Units)');
grid on;