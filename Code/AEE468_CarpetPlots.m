clc; clear; close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%% 1. CONSTANTS
W_payload  = 30800;       % lb  (ISR + munitions + crew)
Wf_W0      = 0.42;        % fuel fraction (patrol-class floor)
C_ew       = -0.07;       % Raymer empty-weight exponent
A_nw       = 0.3977;      % non-wing empty weight coefficient
                           % (back-calculated to hit W0=95,440 at design pt)

% Wing geometry — Raymer Eq. 15.25
AR_wing    = 11.5;
t_c        = 0.12;
lambda     = 0.35;
sweep_qc   = 25 * pi/180;
n_ult      = 3.75;
Scsw_frac  = 0.20;

% Engine weight sensitivity to T/W
TW_ref          = 0.298;
eng_sensitivity = 1.5;    % fractional A_nw increase per unit T/W above baseline

%% 2. SWEEP DEFINITION
WS_fine = linspace(85, 140, 80);   % lb/ft²

% Iso-lines to draw
TW_lines = [0.199, 0.248, 0.298, 0.350];

% Vertical reference lines
WS_refs  = [94, 100, 113, 120, 136];

% Design point
TW_dp = 0.298;
WS_dp = 113;

%% 3. SIZE W0 FOR EACH (T/W, W/S)
nearest = @(v,x) find(abs(v-x)==min(abs(v-x)), 1);

W0_lines = zeros(length(TW_lines), length(WS_fine));

for i = 1:length(TW_lines)
    TW_i = TW_lines(i);
    A_i  = A_nw * (1 + eng_sensitivity*(TW_i - TW_ref));

    for j = 1:length(WS_fine)
        WS_j = WS_fine(j);
        W0   = 95000;   % initial guess

        for iter = 1:3000
            S    = W0 / WS_j;
            Ww   = 0.0051 * (W0*n_ult)^0.557 * S^0.649 * AR_wing^0.5 ...
                   * t_c^(-0.4) * (1+lambda)^0.1 / cos(sweep_qc) ...
                   * (Scsw_frac*S)^0.1;

            denom = 1 - A_i*W0^C_ew - Ww/W0 - Wf_W0;
            if denom < 0.03; W0 = NaN; break; end

            W0_new = W_payload / denom;
            if abs(W0_new - W0)/max(W0,1) < 1e-10
                W0 = W0_new; break;
            end
            W0 = 0.25*W0 + 0.75*W0_new;
        end
        W0_lines(i,j) = W0;
    end
end

% Design point value
di   = nearest(TW_lines, TW_dp);
dj   = nearest(WS_fine,  WS_dp);
W0_dp = W0_lines(di, dj);

fprintf('Design point: W/S = %.0f lb/ft^2 | T/W = %.3f | W0 = %.0f lb (%.1f tons)\n', ...
        WS_dp, TW_dp, W0_dp, W0_dp/2000);

%% 4. PLOT
colors = {[0.00 0.45 0.74],   ...   % blue
          [0.17 0.63 0.17],   ...   % green
          [0.84 0.15 0.15],   ...   % red
          [0.75 0.00 0.75]};        % magenta

figure('Name','2D Carpet Plot','Color','w');
hold on; grid on; box on;

% Iso-T/W lines
for i = 1:length(TW_lines)
    lw = 2.0 + 0.5*(TW_lines(i) == TW_dp);   % design line slightly thicker
    plot(WS_fine, W0_lines(i,:)/1e4, ...
         'Color', colors{i}, 'LineWidth', lw, ...
         'DisplayName', sprintf('T/W = %.3f', TW_lines(i)));
end

% Vertical reference lines
yl = ylim;
for k = 1:length(WS_refs)
    xline(WS_refs(k), '--', 'Color', [0.50 0.50 0.50], ...
          'LineWidth', 1.1, 'HandleVisibility', 'off');
end

% Design point marker
plot(WS_dp, W0_dp/1e4, 'o', ...
     'MarkerSize', 11, ...
     'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', ...
     'LineWidth', 1.8, ...
     'HandleVisibility', 'off');

text(WS_dp + 1.5, W0_dp/1e4 - 0.030, ...
     sprintf('Design Point\nW/S = %.0f lb/ft^2\nT/W = %.3f\nW_0 = %.0f lb', ...
             WS_dp, TW_dp, W0_dp), ...
     'FontSize', 9, 'FontWeight', 'bold', ...
     'BackgroundColor', 'yellow', 'EdgeColor', 'k', 'Margin', 3);

% W/S reference labels along the top
ax = gca;
for k = 1:length(WS_refs)
    text(WS_refs(k) + 0.3, ax.YLim(2)*0.9995, ...
         sprintf('%g', WS_refs(k)), ...
         'FontSize', 8, 'Color', [0.40 0.40 0.40], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

xlabel('W/S  [lb/ft^2]',           'FontSize', 13, 'FontWeight', 'bold');
ylabel('W_0  [lb]  (\times10^4)',   'FontSize', 13, 'FontWeight', 'bold');
title('2D Carpet Plot — ASW Aircraft Conceptual Sizing', ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
xlim([85, 140]);
set(gca, 'FontSize', 11);