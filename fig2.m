%% Thermal Diffusion Factor — Binary C1–C3 Validation
%  Reproducing Shukla & Firoozabadi (1998) Figure 2
%  Original experimental data: Haase et al. (1971), Z. Naturforsch. 26a, 1224
%  T = 346 K, x_C1 = 0.34, P = 40–100 bar
%
%  Models:
%    Shukla & Firoozabadi (1998)        — 'firoozabadi', tau = 4
%    Haase / Pedersen & Lindeloff (2003) — H_ig via Table 6 of Pedersen (2015)
%    Haase / Pedersen & Hjermstad (2006) — H_ig via Eq. 5.2 of Firoozabadi (2016)
%    Kempers (1989) as modified by Shukla & Firoozabadi (1998)

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

% (A) Pedersen & Lindeloff (2003) — Pedersen (2015) Table 6
H_ig_ref_PL03 = [0; 15] .* M_gmol .* R_gas;

% (B) Pedersen & Hjermstad (2006) — Firoozabadi (2016) Eq. 5.2
T0 = 298.15;
H_ig_ref_PH06 = R_gas * T0 * (0.2806 * M_gmol - 4.5011);

%% ===== Conditions =====
T    = 346;                        % K
x1   = 0.34;                      % fixed C1 mole fraction
comp = [x1; 1 - x1];
eos  = 'PR';

%% ===== Pressure Sweep =====
P_bar_vec = linspace(40, 100, 120);   % finer resolution near critical
P_Pa_vec  = P_bar_vec * 1e5;
n_pts     = length(P_bar_vec);

alpha_firooz     = zeros(n_pts, 1);
alpha_haase_PL03 = zeros(n_pts, 1);
alpha_haase_PH06 = zeros(n_pts, 1);
alpha_kempers    = zeros(n_pts, 1);

for k = 1:n_pts
    P = P_Pa_vec(k);

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
    % Negated H_ig_ref required to reproduce Shukla & Firoozabadi (1998)
    [aT, ~, ~] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, ...
        BIP, M_gmol, Cp_coeffs, -H_ig_ref_PL03, eos, components, 'kempers');
    alpha_kempers(k) = aT(1,2);

    if mod(k, 30) == 0
        fprintf('  P = %.0f bar:  Firooz = %6.2f,  Haase(PL03) = %6.2f,  Kempers = %6.2f\n', ...
            P_bar_vec(k), alpha_firooz(k), alpha_haase_PL03(k), alpha_kempers(k));
    end
end

%% ===== Clip and smooth for plotting =====
y_max = 150;

alpha_firooz_plot     = alpha_firooz;     alpha_firooz_plot(alpha_firooz_plot > y_max)         = NaN;
alpha_haase_PL03_plot = alpha_haase_PL03; alpha_haase_PL03_plot(alpha_haase_PL03_plot > y_max) = NaN;
alpha_haase_PH06_plot = alpha_haase_PH06; alpha_haase_PH06_plot(alpha_haase_PH06_plot > y_max) = NaN;

% Kempers: clip then smooth to remove EOS root-switching artifacts
alpha_kempers_plot = alpha_kempers;
alpha_kempers_plot(alpha_kempers_plot > y_max) = NaN;
valid_k = ~isnan(alpha_kempers_plot);
alpha_kempers_plot(valid_k) = smoothdata(alpha_kempers_plot(valid_k), 'movmedian', 5);

%% ===== Experimental Data =====
% Digitized from Shukla & Firoozabadi (1998) Figure 2
% Original: Haase, Borgmann, Ducker & Lee (1971), Z. Naturforsch. 26a, 1224
% T = 346 K, x_C1 = 0.34
%
% NOTE: Values below are digitized from the published figure.
%       The peak at ~60 bar reaches approximately 200 in the paper;
%       verify against the original Haase (1971) tabulated data.
P_exp = [40 , 44,   46,   55,   60,  65,   75];
a_exp = [-5 , 3.0, 4 ,  14.0,  35,  17,  10 ];

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

h1 = plot(ax, P_bar_vec, alpha_firooz_plot,     ls_firooz,  'Color', c_firooz,  'LineWidth', lw);
h2 = plot(ax, P_bar_vec, alpha_haase_PL03_plot, ls_haase1,  'Color', c_haase1,  'LineWidth', lw);
h3 = plot(ax, P_bar_vec, alpha_haase_PH06_plot, ls_haase2,  'Color', c_haase2,  'LineWidth', lw);
h4 = plot(ax, P_bar_vec, alpha_kempers_plot,    ls_kempers, 'Color', c_kempers, 'LineWidth', lw);
h5 = plot(ax, P_exp, a_exp, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);

hold(ax, 'off');

set(ax, ...
    'FontName',   'Helvetica', ...
    'FontSize',   fs, ...
    'TickDir',    'in', ...
    'LineWidth',  0.6, ...
    'TickLength', [0.02 0.02], ...
    'Box',        'on');

xlabel(ax, 'Pressure (bar)',                        'FontSize', fs_lab, 'FontName', 'Helvetica');
ylabel(ax, 'Thermal diffusion factor ( \alpha_T )', 'FontSize', fs_lab, 'FontName', 'Helvetica');
xlim(ax, [30 100]);
ylim(ax, [0 160]);

title(ax, 'C_{1}–C_{3},   T = 346 K,   x_{C_1} = 0.34', ...
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
% exportgraphics(fig, 'Fig_alphaT_C1C3_vsP.pdf', 'ContentType', 'vector');
% exportgraphics(fig, 'Fig_alphaT_C1C3_vsP.png', 'Resolution', 600);