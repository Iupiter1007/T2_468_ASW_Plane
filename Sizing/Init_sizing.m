clc; clear all; close all;
% Aircraft Given Specs
W.crew = 800;   % lb
W.payload = 3e4;    % lb

cruise.C = 0.5;    % specific fuel consumption, lb/hr/lb
loiter.C = 0.4;
S_wet_S_ref = 6;
L_D.max = 17;  % lift / drag
L_D.cruise = 0.866 * L_D.max;

cruise.dist = 2000; % mi, one leg/2
R = cruise.dist/2*5280; % ft
cruise.vel = sqrt(1.4*287*227) * 0.85 * 2.23; % to mph

loiter.time = 4;    % hr
loiter.dist = loiter.time * vel;



WF.fuel = 1.06*(1-W.x/W0)
WF.empty = 0.45 * 0.95;

W2_W0 = 0.97*0.985
W3_W2 = exp(R*cruise.C/cruise.vel/cruise.L_D)

W0 = (W.crew + W.payload) / (1-(WF.fuel)-(WF.empty));

