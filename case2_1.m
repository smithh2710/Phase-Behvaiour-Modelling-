clear; clc; close all;

R = 8.3144598;

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7+'};

comp = [0.1; 2.4; 82.7; 6.43; 2.34; 0.38; 0.65; 0.19; 0.26; 0.35; 4.2] / 100;
M_gmol = [28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; 162.4];
Tc = [126.20; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.40; 628.789];
Pc = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; 20.73] * 1e5;
acentric = [0.04; 0.2250; 0.008; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.296; 0.6542];
ncomp = length(comp);

BIP = zeros(ncomp);
BIP(1,2)= -0.0315; BIP(2,1)= -0.0315;
BIP(1,3)= 0.0278;  BIP(3,1)= 0.0278;
BIP(1,4)= 0.0407;  BIP(4,1)= 0.0407;
BIP(1,5)= 0.0763;  BIP(5,1)= 0.0763;
BIP(1,6)= 0.0944;  BIP(6,1)= 0.0944;
BIP(1,7)= 0.0700;  BIP(7,1)= 0.0700;
BIP(1,8)= 0.0867;  BIP(8,1)= 0.0867;
BIP(1,9)= 0.0878;  BIP(9,1)= 0.0878;
BIP(1,10)= 0.08;   BIP(10,1)= 0.08;
BIP(1,11)= 0.08;   BIP(11,1)= 0.08;
for j = 3:ncomp
    BIP(2,j) = 0.12; BIP(j,2) = 0.12;
end

Zc_all = 0.29056 - 0.08775 * acentric;
Vc_all = Zc_all .* R .* Tc ./ Pc * 1e6;
vt_params.Zc = Zc_all;
vt_params.Vc = Vc_all;
vt_params.components = components;

%% Experimental density data - Regueira et al. (2020), Table 3
%  HPHT Gas Condensate, density in g/cm3
%  NaN = below saturation pressure (no single-phase measurement)

T_exp = [298.15; 323.15; 348.15; 363.15; 373.15; 393.15; 423.75; 463.15];
P_exp_MPa = [60; 80; 100; 120; 140];
nT = length(T_exp);
nP = length(P_exp_MPa);

rho_exp = [
    NaN     0.4476  0.4681  0.4848  0.4987
    NaN     0.4288  0.4510  0.4689  0.4839
    0.3781  0.4107  0.4346  0.4536  0.4694
    0.3649  0.3986  0.4236  0.4434  0.4599
    0.3577  0.3930  0.4189  0.4392  0.4560
    0.3435  0.3809  0.4082  0.4294  0.4466
    0.3226  0.3621  0.3911  0.4135  0.4316
    0.2962  0.3384  0.3693  0.3925  0.4112
];

%% Method definitions

vt_methods = [7, 8, 9];
vt_names = {'SRK', 'SRK + Pina-Martinez', 'SRK + Chen-Li'};
nM = length(vt_methods);

%% Compute model densities at experimental points

rho_model = zeros(nT, nP, nM);
for iT = 1:nT
    for iP = 1:nP
        press = P_exp_MPa(iP) * 1e6;
        Temp  = T_exp(iT);
        for iM = 1:nM
            [rho_kg,~,~,~] = calculate_density(comp, press, Temp, Pc, Tc, acentric, BIP, M_gmol, vt_methods(iM), vt_params);
            rho_model(iT, iP, iM) = rho_kg / 1000;
        end
    end
end

%% AARD (computed only where experimental data exists)

AARD = zeros(nM, 1);
for iM = 1:nM
    valid = ~isnan(rho_exp);
    rel_err = abs(rho_model(:,:,iM) - rho_exp) ./ rho_exp;
    AARD(iM) = mean(rel_err(valid)) * 100;
end

fprintf('%-30s  AARD [%%]\n', 'Method');
fprintf('%s\n', repmat('-', 1, 45));
for iM = 1:nM
    fprintf('%-30s  %.2f\n', vt_names{iM}, AARD(iM));
end

%% Precompute smooth curves

T_plot = [323.15, 348.15, 373.15, 463.15];
T_idx  = [2, 3, 5, 8];
panel_labels = {'(a)', '(b)', '(c)', '(d)'};
nPfine = 60;
P_fine = linspace(55, 145, nPfine);

rho_curves = zeros(4, nPfine, nM);
for k = 1:4
    for ip = 1:nPfine
        for iM = 1:nM
            [rho_kg,~,~,~] = calculate_density(comp, P_fine(ip)*1e6, T_plot(k), Pc, Tc, acentric, BIP, M_gmol, vt_methods(iM), vt_params);
            rho_curves(k, ip, iM) = rho_kg / 1000;
        end
    end
end

%% Style definitions

colors = [
    0.55  0.55  0.55      % SRK no VT    (grey)
    0.00  0.45  0.74      % Pina-Martinez (blue)
    0.47  0.67  0.19      % Chen-Li      (green)
];

line_styles  = {'--', '-.', '-'};
line_widths  = [2.5, 1.8, 1.8];

%% Create figure

fig = figure('Units', 'centimeters', 'Position', [2 2 32 22], ...
             'Color', 'w', 'PaperPositionMode', 'auto');

gap_h = 0.12;
gap_v = 0.10;
marg_l = 0.08;
marg_r = 0.03;
marg_b = 0.14;
marg_t = 0.04;

pw = (1 - marg_l - marg_r - gap_h) / 2;
ph = (1 - marg_b - marg_t - gap_v) / 2;

positions = [
    marg_l,              marg_b + ph + gap_v, pw, ph
    marg_l + pw + gap_h, marg_b + ph + gap_v, pw, ph
    marg_l,              marg_b,              pw, ph
    marg_l + pw + gap_h, marg_b,              pw, ph
];

for k = 1:4
    ax = axes('Position', positions(k,:));
    hold on; box on;

    valid = ~isnan(rho_exp(T_idx(k), :));
    h_exp = plot(P_exp_MPa(valid), rho_exp(T_idx(k), valid), 'ko', ...
        'MarkerSize', 7, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);

    h_lines = gobjects(nM, 1);
    for iM = 1:nM
        h_lines(iM) = plot(P_fine, squeeze(rho_curves(k,:,iM)), ...
            line_styles{iM}, 'Color', colors(iM,:), ...
            'LineWidth', line_widths(iM));
    end

    xlim([50 150]);
    ylim([0.20 0.55]);

    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
            'TickDir', 'in', 'TickLength', [0.015 0.015], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'LineWidth', 0.8);

    if k >= 3
        xlabel('{\itp} / MPa', 'FontSize', 12, 'FontName', 'Times New Roman');
    else
        set(ax, 'XTickLabel', []);
    end

    ylabel('\rho / g\cdotcm^{-3}', 'FontSize', 12, 'FontName', 'Times New Roman');

    text(0.03, 0.95, [panel_labels{k} '  {\itT} = ' sprintf('%.1f K', T_plot(k))], ...
        'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top');

    if k == 1
        h_leg = [h_exp; h_lines];
    end
end

lg = legend(h_leg, [{'Exp. data'}, vt_names], ...
    'Orientation', 'horizontal', 'FontSize', 9, ...
    'FontName', 'Times New Roman', 'NumColumns', 4, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [20, 8];
set(lg, 'Units', 'normalized');
lg_pos = lg.Position;
lg_pos(1) = (1 - lg_pos(3)) / 2;
lg_pos(2) = 0.01;
set(lg, 'Position', lg_pos);