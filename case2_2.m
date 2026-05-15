clear; clc; close all;

R = 8.3144598;

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7+'};
comp = [0.1; 2.4; 82.7; 6.43; 2.34; 0.38; 0.65; 0.19; 0.26; 0.35; 4.2] / 100;
M_gmol = [28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; 162.4];
ncomp = length(comp);

%% PR characterization
Tc_pr = [126.20; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.40; 642.634];
Pc_pr = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; 20.41] * 1e5;
w_pr  = [0.04; 0.2250; 0.008; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.296; 0.5468];

BIP_pr = zeros(ncomp);
BIP_pr(1,2)= -0.0170; BIP_pr(2,1)= -0.0170;
BIP_pr(1,3)= 0.0311;  BIP_pr(3,1)= 0.0311;
BIP_pr(1,4)= 0.0515;  BIP_pr(4,1)= 0.0515;
BIP_pr(1,5)= 0.0852;  BIP_pr(5,1)= 0.0852;
BIP_pr(1,6)= 0.1033;  BIP_pr(6,1)= 0.1033;
BIP_pr(1,7)= 0.0800;  BIP_pr(7,1)= 0.0800;
BIP_pr(1,8)= 0.0922;  BIP_pr(8,1)= 0.0922;
BIP_pr(1,9)= 0.100;   BIP_pr(9,1)= 0.100;
BIP_pr(1,10)= 0.08;   BIP_pr(10,1)= 0.08;
BIP_pr(1,11)= 0.08;   BIP_pr(11,1)= 0.08;
for j = 3:ncomp
    BIP_pr(2,j) = 0.12; BIP_pr(j,2) = 0.12;
end
BIP_pr(2,11) = 0.1000; BIP_pr(11,2) = 0.1000;

% Zc_pr = 0.29056 - 0.08775 * w_pr;
% Vc_pr = Zc_pr .* R .* Tc_pr ./ Pc_pr * 1e6;
Vc_pr = [89.8 ; 94 ; 99 ; 148 ; 203 ; 263 ; 255 ; 306 ; 304 ; 370 ; 690.79]; 
Zc_pr = Pc_pr .* (Vc_pr * 1e-6) ./ (R .* Tc_pr);
ZRA_pr = 0.29056 - 0.08775 * w_pr;
c_peneloux_pr = 0.50033 * R * Tc_pr .* (0.25969 - ZRA_pr) ./ Pc_pr * 1e6;

vtp_pr.Zc = Zc_pr;
vtp_pr.Vc = Vc_pr;
vtp_pr.components = components;
vtp_pr.c_custom = c_peneloux_pr;

%% Method definitions

all_methods = [0, 4, 5, 6];
all_names = {'PR', 'PR + Baled', 'PR + Abudour', 'PR-VT'};
nM = length(all_methods);

%% Experimental data - Regueira et al. (2020), Table 3
%  HPHT Gas Condensate, NaN = below saturation pressure

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

%% Compute model densities at experimental points

rho_model = NaN(nT, nP, nM);
for iT = 1:nT
    for iP = 1:nP
        if isnan(rho_exp(iT, iP)), continue; end
        press = P_exp_MPa(iP) * 1e6;
        Temp  = T_exp(iT);
        for iM = 1:nM
            try
                [rho_kg,~,~,~] = calculate_density(comp, press, Temp, Pc_pr, Tc_pr, w_pr, BIP_pr, M_gmol, all_methods(iM), vtp_pr);
                rho_model(iT, iP, iM) = rho_kg / 1000;
            catch
                rho_model(iT, iP, iM) = NaN;
            end
        end
    end
end

%% AARD (only where experimental data exists)

AARD = zeros(nM, 1);
for iM = 1:nM
    valid = ~isnan(rho_exp) & ~isnan(rho_model(:,:,iM));
    rel_err = abs(rho_model(:,:,iM) - rho_exp) ./ rho_exp;
    AARD(iM) = mean(rel_err(valid)) * 100;
    fprintf('%-30s  AARD = %.2f%%  (%d/%d points valid)\n', ...
        all_names{iM}, AARD(iM), sum(valid(:)), sum(~isnan(rho_exp(:))));
end

%% Precompute smooth curves

T_plot = [348.15, 363.15, 423.75, 463.15];
T_idx  = [3, 4, 7, 8];
panel_labels = {'(a)', '(b)', '(c)', '(d)'};
nPfine = 60;
P_start = [62, 60, 52, 48];

rho_curves = NaN(4, nPfine, nM);
P_fine_all = zeros(4, nPfine);
for k = 1:4
    P_fine_all(k,:) = linspace(P_start(k), 140, nPfine);
    for ip = 1:nPfine
        for iM = 1:nM
            try
                [rho_kg,~,~,~] = calculate_density(comp, P_fine_all(k,ip)*1e6, T_plot(k), Pc_pr, Tc_pr, w_pr, BIP_pr, M_gmol, all_methods(iM), vtp_pr);
                rho_curves(k, ip, iM) = rho_kg / 1000;
            catch
                rho_curves(k, ip, iM) = NaN;
            end
        end
    end
end

%% Style definitions

colors = [
    0.55  0.55  0.55      % PR              (grey)
    0.93  0.69  0.13      % PR+Baled        (gold)
    0.85  0.33  0.10      % PR+Abudour      (orange-red)
    0.49  0.18  0.56      % PR-VT           (purple)
];

line_styles  = {'-', '-', '-', '-'};
line_widths  = [2.5, 1.4, 1.4, 1.4];

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
        h_lines(iM) = plot(P_fine_all(k,:), squeeze(rho_curves(k,:,iM)), ...
            line_styles{iM}, 'Color', colors(iM,:), ...
            'LineWidth', line_widths(iM));
    end

    xlim([52 145]);
    ylim([0.20 0.55]);

    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
            'TickDir', 'in', 'TickLength', [0.015 0.015], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'LineWidth', 0.8);

    xlabel('{\itp} / MPa', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('\rho / g\cdotcm^{-3}', 'FontSize', 12, 'FontName', 'Times New Roman');

    text(0.03, 0.95, [panel_labels{k} '  {\itT} = ' sprintf('%.1f K', T_plot(k))], ...
        'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top');

    if k == 1
        h_leg = [h_exp; h_lines];
    end
end

lg = legend(h_leg, [{'Exp. data'}, all_names], ...
    'Orientation', 'horizontal', 'FontSize', 8.5, ...
    'FontName', 'Times New Roman', 'NumColumns', 5, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [20, 8];
set(lg, 'Units', 'normalized');
lg_pos = lg.Position;
lg_pos(1) = (1 - lg_pos(3)) / 2;
lg_pos(2) = 0.005;
set(lg, 'Position', lg_pos);