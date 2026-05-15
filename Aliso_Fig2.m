clc; clear; close all;

R = 8.3144598;

components = {'H2', 'CH4', 'C2H6'};
n = 3;

comp_ref = [0.10; 0.85; 0.05];

M_gmol    = [2.016;  16.043;  30.070];
Pc        = [12.90;  46.00;   48.72 ] * 1e5;
Tc        = [33.15;  190.60;  305.32];
acentric  = [-0.22;  0.011;   0.0995];

Zc = 0.2905 - 0.085 * acentric;

Cp_coeffs = [
    29.11,  -1.916e-3,   4.003e-6,  -8.704e-10;
    19.25,   5.213e-2,   1.197e-5,  -1.132e-8;
     5.409,  1.781e-1,  -6.938e-5,   8.713e-9;
];

H_ig_ref = [0; 0; 0];

BIP = zeros(n, n);

h_ref     = 2351.3;
h_top     = h_ref - 61;
h_bot     = h_ref + 61;
press_ref = 310.26e5;
temp_ref  = 352.32;

tol     = 1e-6;
maxiter = 500;

vt_method = 9;
vt_params = struct();
vt_params.Zc         = Zc;
vt_params.components = components;

dTdh_vals   = [0, 0.010, 0.020, 0.027, 0.040, 0.050];
dTdh_labels = {'Isothermal', 'dT/dh = 0.010 K/m', 'dT/dh = 0.020 K/m', ...
               'dT/dh = 0.027 K/m (Aliso Canyon)', 'dT/dh = 0.040 K/m', 'dT/dh = 0.050 K/m'};
ng = length(dTdh_vals);

line_styles = {'-', '--', ':', '-.', '--', '-'};
line_widths = [2.0, 2.0, 2.0, 2.6, 2.0, 2.0];
line_colors = {
    [0.00 0.45 0.74],   % blue       - isothermal
    [0.00 0.60 0.50],   % teal       - 0.010
    [0.93 0.69 0.13],   % amber      - 0.020
    [0.85 0.10 0.10],   % red        - 0.027 Aliso Canyon
    [0.93 0.55 0.13],   % orange     - 0.040
    [0.58 0.30 0.71]    % purple     - 0.050
};

h_vec = linspace(h_top, h_bot, 100)';
nd    = length(h_vec);

xH2_all = zeros(nd, ng);
P_all   = zeros(nd, ng);
rho_all = zeros(nd, ng);

for g = 1:ng
    dTdh = dTdh_vals(g);
    for i = 1:nd
        [comp_h, P_h, T_h, ~, ~] = main_hasse(h_vec(i), h_ref, comp_ref, press_ref, temp_ref, ...
            dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        xH2_all(i,g) = comp_h(1) * 100;
        P_all(i,g)   = P_h / 1e5;
        try
            rho_all(i,g) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
        catch
            rho_all(i,g) = NaN;
        end
    end
end

fn       = 'Times New Roman';
fs_tick  = 13;
fs_lab   = 14;
fs_pan   = 15;
fs_leg   = 12;

pw    = 9;
ph    = 13;
gap_h = 1.5;
ml = 2.0; mr = 0.6; mt = 0.6; mb = 4.5;

fig = figure;
fig.Units    = 'centimeters';
fig.Color    = 'w';
fig.Position = [2, 2, 3*pw + 2*gap_h + ml + mr, ph + mt + mb];
fig.PaperPositionMode = 'auto';

x_data    = {xH2_all,  rho_all,           P_all};
xlabels   = {'H_{2} (mol%)', 'Density (kg/m^{3})', 'Pressure (bar)'};
panlabels = {'(a)', '(b)', '(c)'};

for k = 1:3
    left_cm = ml + (k-1)*(pw + gap_h);
    ax = axes('Units', 'centimeters', 'Position', [left_cm, mb, pw, ph]);
    hold on;
    h_lines = gobjects(ng, 1);
    for g = 1:ng
        h_lines(g) = plot(x_data{k}(:,g), h_vec, ...
            'LineStyle', line_styles{g}, ...
            'Color',     line_colors{g}, ...
            'LineWidth', line_widths(g));
    end
    yline(h_ref, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    hold off;
    set(ax, 'YDir', 'reverse', 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'FontName', fn, 'FontSize', fs_tick, 'Layer', 'top', ...
        'LineWidth', 0.9, 'TickLength', [0.015 0.015]);
    xlabel(xlabels{k}, 'FontName', fn, 'FontSize', fs_lab);
    ylabel('Depth (m)', 'FontName', fn, 'FontSize', fs_lab);
    text(0.05, 0.97, panlabels{k}, 'Units', 'normalized', ...
        'FontName', fn, 'FontSize', fs_pan, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    if k == 2
        leg = legend(h_lines, dTdh_labels, ...
            'Orientation', 'horizontal', ...
            'NumColumns', 3, ...
            'FontName', fn, 'FontSize', fs_leg, ...
            'Box', 'on', 'EdgeColor', [0.4 0.4 0.4], ...
            'LineWidth', 0.7, 'Color', 'w');
        leg.ItemTokenSize = [28, 14];
        leg.Units    = 'centimeters';
        leg.Position(1) = ml;
        leg.Position(2) = 0.4;
        leg.Position(3) = 3*pw + 2*gap_h;
    end
end