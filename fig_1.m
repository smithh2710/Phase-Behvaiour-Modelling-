%% Thermal Diffusion Factor — Binary C1–C3 Validation
%  Reproducing Shukla & Firoozabadi (1998) Figure 3
%  Original experimental data: Haase et al. (1971), Z. Naturforsch. 26a, 1224
%  T = 346 K, P = 55 bar, x_C1 = 0.30–0.70

clear; clc; close all;

%% ===== Component Properties =====
components = {'C1', 'C3'};
M_gmol   = [16.043; 44.096];
Tc       = [190.56; 369.83];
Pc       = [45.99; 42.48] * 1e5;
acentric = [0.0115; 0.1523];
BIP      = zeros(2);

Cp_coeffs = [
    19.25,   0.05213,   1.197e-5,  -1.132e-8;
    -4.224,  0.3063,   -1.586e-4,   3.215e-8];

%% ===== Ideal-Gas Reference Enthalpy Parameterizations =====
R_gas = 8.3144598;

H_ig_ref_PL03 = [0; 15] .* M_gmol .* R_gas;

T0 = 298.15;
H_ig_ref_PH06 = R_gas * T0 * (0.2806 * M_gmol - 4.5011);

%% ===== Conditions =====
T = 346;
P = 55e5;
eos = 'PR';

%% ===== Composition Sweep =====
x_C1_vec = linspace(0.30, 0.70, 50);
n_pts    = length(x_C1_vec);

alpha_firooz     = zeros(n_pts, 1);
alpha_haase_PL03 = zeros(n_pts, 1);
alpha_haase_PH06 = zeros(n_pts, 1);
alpha_kempers    = zeros(n_pts, 1);

for k = 1:n_pts
    x1   = x_C1_vec(k);
    comp = [x1; 1 - x1];

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'firoozabadi', ...
        'tau', 4.0);
    alpha_firooz(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'haase');
    alpha_haase_PL03(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PH06, eos, components, 'haase');
    alpha_haase_PH06(k) = aT(1,2);

    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, -H_ig_ref_PL03, eos, components, 'kempers');
    alpha_kempers(k) = aT(1,2);
end

%% ===== Clip Kempers for plotting =====
alpha_kempers_plot = alpha_kempers;
alpha_kempers_plot(alpha_kempers_plot > 23) = NaN;

%% ===== Experimental Data =====
x_exp = [0.34,  0.35,  0.42,  0.51,  0.58,  0.63];
a_exp = [12.0,  12.0,   2.7,   1.0,   0.6,   0.8];

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

h1 = plot(ax, x_C1_vec, alpha_firooz,      '-',  'Color', c_firooz,  'LineWidth', lw);
h2 = plot(ax, x_C1_vec, alpha_haase_PL03,  '--', 'Color', c_haase1,  'LineWidth', lw);
h3 = plot(ax, x_C1_vec, alpha_haase_PH06,  '-.', 'Color', c_haase2,  'LineWidth', lw);
h4 = plot(ax, x_C1_vec, alpha_kempers_plot,':',  'Color', c_kempers, 'LineWidth', lw);
h5 = plot(ax, x_exp, a_exp, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);

hold(ax, 'off');

set(ax, ...
    'FontName',   'Times New Roman', ...
    'FontSize',   11, ...
    'TickDir',    'in', ...
    'TickLength', [0.015 0.015], ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'LineWidth',  0.8);

xlabel(ax, 'Mole fraction of methane', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax, 'Thermal diffusion factor ({\it\alpha}_{T})', 'FontSize', 12, 'FontName', 'Times New Roman');
xlim(ax, [0.28 0.72]);
ylim(ax, [-5 25]);

lgd = legend(ax, [h5, h1, h2, h3, h4], ...
    {'Exp. data (Haase et al., 1971)', ...
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
% exportgraphics(fig, 'Fig_alphaT_C1C3_55bar.pdf', 'ContentType', 'vector');
% exportgraphics(fig, 'Fig_alphaT_C1C3_55bar.png', 'Resolution', 600);