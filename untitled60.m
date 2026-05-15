%% Thermal Diffusion Factor — Binary C7–C16 Validation
%  Reproducing Shukla & Firoozabadi (1998) Figure 12
%  Original experimental data: Shieh (1969), J. Phys. Chem. 73, 1508
%  T = 308 K, P = 1.01325 bar, x_C7 = 0.0–1.0
%
%  Models:
%    Shukla & Firoozabadi (1998)        — 'firoozabadi', tau = 4
%    Haase / Pedersen & Lindeloff (2003) — H_ig via Table 6 of Pedersen (2015)
%    Haase / Pedersen & Hjermstad (2006) — H_ig via Eq. 5.2 of Firoozabadi (2016)
%    Kempers (1989) as modified by Shukla & Firoozabadi (1998)

clear; clc; close all;

%% ===== Component Properties =====
components = {'C7', 'C16'};
M_gmol   = [100.204; 226.441];
Tc       = [540.20;  723.00];          % K
Pc       = [27.40;   14.00] * 1e5;    % Pa
acentric = [0.3495;  0.7174];
BIP      = [0, 0.01; 0.01, 0];        % kij = 0.01

% Ideal-gas Cp coefficients (Reid et al.)
Cp_coeffs = [
    -5.146,  0.6762,  -3.651e-4,  7.658e-8;    % C7
    -7.046,  1.4850,  -8.528e-4,  1.939e-7];    % C16

%% ===== Ideal-Gas Reference Enthalpy Parameterizations =====
R_gas = 8.3144598;

% (A) Pedersen & Lindeloff (2003) — Pedersen (2015) Table 6
%     H_ig / (M*R) [K]:  C7 = 105,  C16 = 260  (approximate)
H_ig_ref_PL03 = [105; 260] .* M_gmol .* R_gas;

% (B) Pedersen & Hjermstad (2006) — Firoozabadi (2016) Eq. 5.2
T0 = 298.15;
H_ig_ref_PH06 = R_gas * T0 * (0.2806 * M_gmol - 4.5011);

%% ===== Conditions =====
T = 308;                           % K  (35°C)
P = 1.01325e5;                     % Pa (1 atm)
eos = 'PR';

%% ===== Composition Sweep =====
x_C7_vec = linspace(0.05, 0.95, 80);
n_pts    = length(x_C7_vec);

alpha_firooz     = zeros(n_pts, 1);
alpha_haase_PL03 = zeros(n_pts, 1);
alpha_haase_PH06 = zeros(n_pts, 1);
alpha_kempers    = zeros(n_pts, 1);

for k = 1:n_pts
    x1   = x_C7_vec(k);
    comp = [x1; 1 - x1];

    % --- Shukla & Firoozabadi (1998), tau = 4 ---
    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'firoozabadi', ...
        'tau', 4.0);
    alpha_firooz(k) = aT(1,2);

    % --- Haase: Pedersen & Lindeloff (2003) parameterization ---
    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'haase');
    alpha_haase_PL03(k) = aT(1,2);

    % --- Haase: Pedersen & Hjermstad (2006) parameterization ---
    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PH06, eos, components, 'haase');
    alpha_haase_PH06(k) = aT(1,2);

    % --- Kempers (1989), modified ---
    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, H_ig_ref_PL03, eos, components, 'kempers');
    alpha_kempers(k) = aT(1,2);

    if mod(k, 20) == 0
        fprintf('  x_C7 = %.2f:  Firooz = %6.3f,  Haase(PL03) = %6.3f,  Kempers = %6.3f\n', ...
            x1, alpha_firooz(k), alpha_haase_PL03(k), alpha_kempers(k));
    end
end

%% ===== Experimental Data =====
% Digitized from Shukla & Firoozabadi (1998) Figure 12
% Original: Shieh (1969), J. Phys. Chem. 73, 1508
% T = 308 K (35°C), P = 1 atm
%
% NOTE: alpha values are small (~0.5–2.0) compared to light HC systems.
%       The paper notes that Rutherford (Haase) model predicts wrong sign.
%       Firoozabadi model captures correct sign and magnitude.
x_exp = [0.10,  0.20,  0.30,  0.40,  0.50,  0.60,  0.70,  0.80,  0.90];
a_exp = [0.70,  0.80,  0.85,  0.75,  0.60,  0.50,  0.45,  0.55,  0.80];

%% ===== Figure =====
% Colour palette — consistent with Reservoir 1 figures
c_firooz  = [0.49 0.18 0.56];   % purple
c_haase1  = [0.85 0.33 0.10];   % red
c_haase2  = [0.00 0.45 0.74];   % blue
c_kempers = [0.47 0.67 0.19];   % green

ls_firooz  = '-';
ls_haase1  = '--';
ls_haase2  = '-.';
ls_kempers = ':';

lw     = 1.6;
fs     = 9;
fs_lab = 10;

fig = figure('Units', 'centimeters', 'Position', [3 3 14 10], ...
    'Color', 'w', 'PaperPositionMode', 'auto');
ax = axes(fig);
hold(ax, 'on');

h1 = plot(ax, x_C7_vec, alpha_firooz,     ls_firooz,  'Color', c_firooz,  'LineWidth', lw);
h2 = plot(ax, x_C7_vec, alpha_haase_PL03, ls_haase1,  'Color', c_haase1,  'LineWidth', lw);
h3 = plot(ax, x_C7_vec, alpha_haase_PH06, ls_haase2,  'Color', c_haase2,  'LineWidth', lw);
h4 = plot(ax, x_C7_vec, alpha_kempers,    ls_kempers, 'Color', c_kempers, 'LineWidth', lw);
h5 = plot(ax, x_exp, a_exp, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);

hold(ax, 'off');

set(ax, ...
    'FontName',   'Helvetica', ...
    'FontSize',   fs, ...
    'TickDir',    'in', ...
    'LineWidth',  0.6, ...
    'TickLength', [0.02 0.02], ...
    'Box',        'on');

xlabel(ax, 'C_{7} mole fraction',                   'FontSize', fs_lab, 'FontName', 'Helvetica');
ylabel(ax, 'Thermal diffusion factor ( \alpha_T )',  'FontSize', fs_lab, 'FontName', 'Helvetica');
% xlim(ax, [0.0 1.0]);
% ylim(ax, [-3 3]);

title(ax, 'C_{7}–C_{16},   T = 308 K,   P = 1 atm', ...
    'FontSize', 10, 'FontName', 'Helvetica', 'FontWeight', 'normal');

lgd = legend(ax, [h1, h2, h3, h4, h5], ...
    {'Firoozabadi (1998), \tau = 4', ...
     'Haase — Pedersen & Lindeloff (2003)', ...
     'Haase — Pedersen & Hjermstad (2006)', ...
     'Kempers (1989)', ...
     'Experimental'}, ...
    'FontSize',  8, ...
    'FontName',  'Helvetica', ...
    'Location',  'northeast', ...
    'Box',       'on', ...
    'EdgeColor', [0.5 0.5 0.5]);
lgd.ItemTokenSize = [20, 8];

%% ===== Export =====
% exportgraphics(fig, 'Fig_alphaT_C7C16_308K.pdf', 'ContentType', 'vector');
% exportgraphics(fig, 'Fig_alphaT_C7C16_308K.png', 'Resolution', 600);