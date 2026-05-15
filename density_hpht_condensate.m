clear; clc; close all;

R = 8.3144598;

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7+'};
comp = [0.1; 2.4; 82.7; 6.43; 2.34; 0.38; 0.65; 0.19; 0.26; 0.35; 4.2] / 100;
M_gmol = [28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; 162.4];
ncomp = length(comp);

%% SRK characterization
Tc_srk = [126.20; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.40; 628.789];
Pc_srk = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; 20.73] * 1e5;
w_srk  = [0.04; 0.2250; 0.008; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.296; 0.6542];

BIP_srk = zeros(ncomp);
BIP_srk(1,2)= -0.0315; BIP_srk(2,1)= -0.0315;
BIP_srk(1,3)= 0.0278;  BIP_srk(3,1)= 0.0278;
BIP_srk(1,4)= 0.0407;  BIP_srk(4,1)= 0.0407;
BIP_srk(1,5)= 0.0763;  BIP_srk(5,1)= 0.0763;
BIP_srk(1,6)= 0.0944;  BIP_srk(6,1)= 0.0944;
BIP_srk(1,7)= 0.0700;  BIP_srk(7,1)= 0.0700;
BIP_srk(1,8)= 0.0867;  BIP_srk(8,1)= 0.0867;
BIP_srk(1,9)= 0.0878;  BIP_srk(9,1)= 0.0878;
BIP_srk(1,10)= 0.08;   BIP_srk(10,1)= 0.08;
BIP_srk(1,11)= 0.08;   BIP_srk(11,1)= 0.08;
for j = 3:ncomp
    BIP_srk(2,j) = 0.12; BIP_srk(j,2) = 0.12;
end

Zc_srk = 0.29056 - 0.08775 * w_srk;
vtp_srk.Zc = Zc_srk;
vtp_srk.Vc = Zc_srk .* R .* Tc_srk ./ Pc_srk * 1e6;
vtp_srk.components = components;

%% PR characterization
Tc_pr = [126.20; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.40; 628.789];
Pc_pr = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; 20.73] * 1e5;
w_pr  = [0.04; 0.2250; 0.008; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.296; 0.6542];

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

Zc_pr = 0.29056 - 0.08775 * w_pr;
Vc_pr = Zc_pr .* R .* Tc_pr ./ Pc_pr * 1e6;
ZRA_pr = Zc_pr;
c_peneloux_pr = 0.50033 * R * Tc_pr .* (0.25969 - ZRA_pr) ./ Pc_pr * 1e6;

vtp_pr.Zc = Zc_pr;
vtp_pr.Vc = Vc_pr;
vtp_pr.components = components;
vtp_pr.c_custom = c_peneloux_pr;

%% Method definitions

all_methods = [0, 3, 4, 5, 6, 7, 8, 9];
all_names = {'PR', 'PR + Ungerer-Batut', ...
             'PR + Baled', 'PR + Abudour', 'PR + Constant VT', ...
             'SRK', 'SRK + Pina-Martinez', 'SRK + Chen-Li'};
nM = length(all_methods);
is_srk = all_methods >= 7;

%% Experimental data - Regueira et al. (2020), Table 3

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
                if is_srk(iM)
                    [rho_kg,~,~,~] = calculate_density(comp, press, Temp, Pc_srk, Tc_srk, w_srk, BIP_srk, M_gmol, all_methods(iM), vtp_srk);
                else
                    [rho_kg,~,~,~] = calculate_density(comp, press, Temp, Pc_pr, Tc_pr, w_pr, BIP_pr, M_gmol, all_methods(iM), vtp_pr);
                end
                rho_model(iT, iP, iM) = rho_kg / 1000;
            catch
                rho_model(iT, iP, iM) = NaN;
            end
        end
    end
end

%% AARD per isotherm

AARD_iso   = NaN(nM, nT);
maxErr_iso = NaN(nM, nT);
N_iso      = zeros(1, nT);

for iT = 1:nT
    N_iso(iT) = sum(~isnan(rho_exp(iT,:)));
end

for iM = 1:nM
    for iT = 1:nT
        valid = ~isnan(rho_exp(iT,:)) & ~isnan(rho_model(iT,:,iM));
        if any(valid)
            rel_err = abs(rho_model(iT,valid,iM) - rho_exp(iT,valid)) ./ rho_exp(iT,valid);
            AARD_iso(iM, iT)   = mean(rel_err) * 100;
            maxErr_iso(iM, iT) = max(rel_err)  * 100;
        end
    end
end

AARD_overall   = zeros(nM, 1);
maxErr_overall = zeros(nM, 1);
for iM = 1:nM
    valid = ~isnan(rho_exp) & ~isnan(rho_model(:,:,iM));
    rel_err_all = abs(rho_model(:,:,iM) - rho_exp) ./ rho_exp;
    AARD_overall(iM)   = mean(rel_err_all(valid)) * 100;
    maxErr_overall(iM) = max(rel_err_all(valid)) * 100;
end

fprintf('\n=== AARD per isotherm (%%) ===\n');
fprintf('%-25s', 'Method');
for iT = 1:nT, fprintf(' %8.2fK', T_exp(iT)); end
fprintf(' %10s\n', 'Overall');
fprintf('%-25s', 'N points');
for iT = 1:nT, fprintf(' %9d', N_iso(iT)); end
fprintf(' %10d\n', sum(N_iso));
fprintf('%s\n', repmat('-', 1, 25 + 9*nT + 11));
for iM = 1:nM
    fprintf('%-25s', all_names{iM});
    for iT = 1:nT
        if isnan(AARD_iso(iM, iT))
            fprintf(' %9s', '--');
        else
            fprintf(' %9.2f', AARD_iso(iM, iT));
        end
    end
    fprintf(' %10.2f\n', AARD_overall(iM));
end

fprintf('\n=== Max relative error per isotherm (%%) ===\n');
fprintf('%-25s', 'Method');
for iT = 1:nT, fprintf(' %8.2fK', T_exp(iT)); end
fprintf(' %10s\n', 'Overall');
fprintf('%s\n', repmat('-', 1, 25 + 9*nT + 11));
for iM = 1:nM
    fprintf('%-25s', all_names{iM});
    for iT = 1:nT
        if isnan(maxErr_iso(iM, iT))
            fprintf(' %9s', '--');
        else
            fprintf(' %9.2f', maxErr_iso(iM, iT));
        end
    end
    fprintf(' %10.2f\n', maxErr_overall(iM));
end

fprintf('\nN total: %d valid points\n', sum(N_iso));

%% Predicted density per (T, P) point - for appendix table

col_order  = [1, 6, 2, 3, 5, 4, 7, 8];
col_labels = {'PR','SRK','PR+UB','PR+Baled','PR+Const','PR+Abudour','SRK+P-M','SRK+Chen-Li'};

fprintf('\n=== Predicted density per (T, P) point (g/cm3) ===\n');
fprintf('%-9s %-9s %-9s', 'T (K)', 'P (MPa)', 'rho_exp');
for iC = 1:length(col_order)
    fprintf(' %-12s', col_labels{iC});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 9*3 + 13*length(col_order)));

for iT = 1:nT
    for iP = 1:nP
        if isnan(rho_exp(iT, iP)), continue; end
        fprintf('%-9.2f %-9.0f %-9.4f', T_exp(iT), P_exp_MPa(iP), rho_exp(iT, iP));
        for iC = 1:length(col_order)
            iM = col_order(iC);
            if isnan(rho_model(iT, iP, iM))
                fprintf(' %-12s', '--');
            else
                fprintf(' %-12.4f', rho_model(iT, iP, iM));
            end
        end
        fprintf('\n');
    end
end

%% Precompute smooth curves

T_plot = [348.15, 363.15, 423.75, 463.15];
T_idx  = [3, 4, 7, 8];
nPfine = 60;
P_start = [62, 60, 52, 48];

rho_curves = NaN(4, nPfine, nM);
P_fine_all = zeros(4, nPfine);
for k = 1:4
    P_fine_all(k,:) = linspace(P_start(k), 140, nPfine);
    for ip = 1:nPfine
        for iM = 1:nM
            try
                if is_srk(iM)
                    [rho_kg,~,~,~] = calculate_density(comp, P_fine_all(k,ip)*1e6, T_plot(k), Pc_srk, Tc_srk, w_srk, BIP_srk, M_gmol, all_methods(iM), vtp_srk);
                else
                    [rho_kg,~,~,~] = calculate_density(comp, P_fine_all(k,ip)*1e6, T_plot(k), Pc_pr, Tc_pr, w_pr, BIP_pr, M_gmol, all_methods(iM), vtp_pr);
                end
                rho_curves(k, ip, iM) = rho_kg / 1000;
            catch
                rho_curves(k, ip, iM) = NaN;
            end
        end
    end
end

%% Style definitions

colors = [
    0.55  0.55  0.55
    0.47  0.67  0.19
    0.93  0.69  0.13
    0.85  0.33  0.10
    0.49  0.18  0.56
    0.55  0.55  0.55
    0.30  0.75  0.93
    0.00  0.50  0.00
];

line_styles  = {'-', '-', '-', '-', '-', '--', '--', '--'};
line_widths  = [2.5, 1.4, 1.4, 1.4, 1.4, 2.5, 1.4, 1.4];

%% Create 4 separate figures

for k = 1:4
    fig = figure('Units', 'centimeters', 'Position', [2+2*k 2+1*k 14 10], ...
                 'Color', 'w', 'PaperPositionMode', 'auto');

    ax = axes('Position', [0.13 0.15 0.83 0.80]);
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

    xlim([45 145]);
    ylim([0.20 0.55]);

    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
            'TickDir', 'in', 'TickLength', [0.015 0.015], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'LineWidth', 0.8);

    xlabel('{\itp} / MPa', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('\rho / g\cdotcm^{-3}', 'FontSize', 12, 'FontName', 'Times New Roman');

    text(0.03, 0.97, sprintf('{\\itT} = %.2f K', T_plot(k)), ...
        'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top');

    loc = 'southeast';
    if k >= 3, loc = 'northwest'; end
    lg = legend([h_exp; h_lines], [{'Exp. data'}, all_names], ...
        'Location', loc, 'FontSize', 8, ...
        'FontName', 'Times New Roman', ...
        'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
    lg.ItemTokenSize = [18, 6];
end