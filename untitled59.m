%% Thermal Diffusion Factor — Binary C1–nC4 Validation
%  Reproducing Shukla & Firoozabadi (1998) Figure 6
%  Original experimental data: Rutherford & Roof (1959), J. Phys. Chem. 63, 1506
%  T = 344 K, x_C1 = 0.49, P = 100–207 bar

clear; clc; close all;

%% ===== Component Properties =====
components = {'C1', 'nC4'};
M_gmol   = [16.043; 58.123];
Tc       = [190.56; 425.12];
Pc       = [45.99; 37.96] * 1e5;
acentric = [0.0115; 0.2002];
BIP      = [0, 0.01; 0.01, 0];

Cp_coeffs = [
    19.25,   0.05213,   1.197e-5,  -1.132e-8;
     9.487,  0.3313,   -1.108e-4,  -2.822e-9];

%% ===== Ideal-Gas Reference Enthalpy Parameterizations =====
R_gas = 8.3144598;

H_ig_ref_PL03 = [0; 45] .* M_gmol .* R_gas;

T0 = 298.15;
H_ig_ref_PH06 = R_gas * T0 * (0.2806 * M_gmol - 4.5011);

%% ===== Conditions =====
T    = 344;
x1   = 0.49;
comp = [x1; 1 - x1];
eos  = 'PR';

%% ===== Pressure Sweep =====
P_bar_vec = linspace(100, 210, 100);
P_Pa_vec  = P_bar_vec * 1e5;
n_pts     = length(P_bar_vec);

alpha_firooz     = zeros(n_pts, 1);
alpha_haase_PL03 = zeros(n_pts, 1);
alpha_haase_PH06 = zeros(n_pts, 1);
alpha_kempers    = zeros(n_pts, 1);

for k = 1:n_pts
    P = P_Pa_vec(k);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'firoozabadi', ...
        'tau',4.0);
    alpha_firooz(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'haase');
    alpha_haase_PL03(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PH06, eos, components, 'haase');
    alpha_haase_PH06(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, [0;0], eos, components, 'kempers');
    alpha_kempers(k) = aT(1,2);

    if mod(k, 25) == 0
        fprintf('  P = %.0f bar:  Firooz = %6.2f,  Haase(PL03) = %6.2f,  Kempers = %6.2f\n', ...
            P_bar_vec(k), alpha_firooz(k), alpha_haase_PL03(k), alpha_kempers(k));
    end
end

%% ===== Experimental Data =====
P_exp = [  122,   126,   134,   140,   165,   198];
a_exp = [  8.0,   7,   5.5,   5.0,   4.5,  4 ];

%% ===== Figure =====

c_firooz  = [0.49 0.18 0.56];
c_haase1  = [0.85 0.33 0.10];
c_haase2  = [0.00 0.45 0.74];
c_kempers = [0.00 0.60 0.30];

lw = 2.0;

fig = figure('Units', 'centimeters', 'Position', [3 3 14 10], ...
    'Color', 'w', 'PaperPositionMode', 'auto');
ax = axes(fig);
hold(ax, 'on'); box(ax, 'on');

h1 = plot(ax, P_bar_vec, alpha_firooz,     '-',  'Color', c_firooz,  'LineWidth', lw);
h2 = plot(ax, P_bar_vec, alpha_haase_PL03, '--', 'Color', c_haase1,  'LineWidth', lw);
h3 = plot(ax, P_bar_vec, alpha_haase_PH06, '-.', 'Color', c_haase2,  'LineWidth', lw);
h4 = plot(ax, P_bar_vec, alpha_kempers,    ':',  'Color', c_kempers, 'LineWidth', lw);
h5 = plot(ax, P_exp, a_exp, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);

hold(ax, 'off');

set(ax, ...
    'FontName',   'Times New Roman', ...
    'FontSize',   11, ...
    'TickDir',    'in', ...
    'TickLength', [0.015 0.015], ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'LineWidth',  0.8);

xlabel(ax, 'Pressure (bar)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax, 'Thermal diffusion factor ({\it\alpha}_{T})', 'FontSize', 12, 'FontName', 'Times New Roman');
xlim(ax, [95 210]);
ylim(ax, [-10 40]);

lgd = legend(ax, [h5, h1, h2, h3, h4], ...
    {'Exp. data (Rutherford & Roof, 1959)', ...
     'Firoozabadi (1998), \tau = 4', ...
     'Haase — Pedersen & Lindeloff (2003)', ...
     'Haase — Pedersen & Hjermstad (2006)', ...
     'Kempers (1989)'}, ...
    'FontSize',  9, ...
    'FontName',  'Times New Roman', ...
    'Location',  'northeast', ...
    'Box',       'on', ...
    'EdgeColor', [0.5 0.5 0.5], ...
    'Interpreter', 'tex');
lgd.ItemTokenSize = [22, 10];

%% ===== Export =====
% exportgraphics(fig, 'Fig_alphaT_C1nC4_344K.pdf', 'ContentType', 'vector');
% exportgraphics(fig, 'Fig_alphaT_C1nC4_344K.png', 'Resolution', 600);