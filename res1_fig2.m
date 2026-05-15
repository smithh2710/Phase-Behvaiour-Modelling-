clear; clc; close all;

components = {'N2', 'CO2', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', 'nC6', 'nC7', 'nC8', 'nC9', 'nC10', 'nC12', 'nC14','nC17', 'nC19', 'nC22', 'nC26', 'C30-C37', 'C38-C80'};

comp_ref = [0.42; 0.71; 49.885; 7.771; 6.741; 1.030; 3.21; 1.16; 1.55; 1.88 ;3.490; 3.730 ; 2.31 ; 3.293; 2.620; 2.961; 1.480; 1.673; 1.187; 1.264;0.982; 0.652];
comp_ref = comp_ref ./ sum(comp_ref) ;
M_gmol = [28.01,44.01 ,16.04 ,30.07 ,44.10 ,58.12 ,58.12 ,72.15 ,72.15 ,86.18 ,96 ,107 ,121 ,140.09 ,167.57 ,204.75 ,243.6 ,275.27 ,317.02 ,370.39 ,456.83 ,640.76]';
acentric = [0.0400 ,0.2250 ,0.0080 ,0.0980 ,0.1520 ,0.1760 ,0.1930 ,0.2270 ,0.2510 ,0.2960 ,0.3375 ,0.3743 ,0.4204 ,0.4833 ,0.5703 ,0.6847 ,0.7947 ,0.8805 ,0.9836 ,1.0983 ,1.2292 ,1.1596]' ;
Tc = [-146.95 ,31.05 ,-82.55 ,32.25 ,96.65 ,134.95 ,152.05 ,187.25 ,196.45 234.25 ,301.63 ,323.97 ,349.62 ,381.84 ,422.96 ,473.01 ,519.43 ,555.64 ,600.26 ,654.86 ,737.85 ,912.38]' + 273.15 ;
Pc = [33.94 ,73.76 ,46.00 ,48.84 ,42.46 ,36.48 ,38 ,33.84 ,33.74 ,29.69 ,29.59 ,27.37 ,25.01 ,22.59 ,20.12 ,17.91 ,16.4 ,15.55 ,14.71 ,13.96 ,13.16 ,12.27 ]'  * 1e5 ;
h_ref = 204;
press_ref = 286e5;
temp_ref = 94 + 273.15;
c_JY = [-4.23; -1.64; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12;1.39; 15.53 ; 19.96 ; 25.14 ; 31.89 ; 39.65 ; 46.52 ; 49.62 ; 50.37 ; 48.37; 42.77 ; 27.92; -6.04];

R = 8.3144598;
n = length(comp_ref);
BIP = zeros(n, n);
BIP(1,2) = -0.0170; BIP(2,1) = -0.0170;
BIP(1,3) = 0.0311;  BIP(3,1) = 0.0311;
BIP(1,4) = 0.0515;  BIP(4,1) = 0.0515;
BIP(1,5) = 0.0852;  BIP(5,1) = 0.0852;
BIP(1,6) = 0.1033;  BIP(6,1) = 0.1033;
BIP(1,7) = 0.0800;  BIP(7,1) = 0.0800;
BIP(1,8) = 0.0922;  BIP(8,1) = 0.0922;
BIP(1,9) = 0.1000;  BIP(9,1) = 0.1000;
for i = 10:n, BIP(1,i) = 0.0800; BIP(i,1) = 0.0800; end
for i = 3:10, BIP(2,i) = 0.1200; BIP(i,2) = 0.1200; end
for i = 11:n, BIP(2,i) = 0.100;  BIP(i,2) = 0.100;  end

Vc = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0; 432; 492; 548; 600.22; 697.08; 810.90; 910.91; 980.42; 1061.51;1152.43; 1287.26; 1648.29];
 % Zc = 0.29056 - 0.08775 .* acentric ;
   Zc = Pc .* (Vc * 1e-6) ./ (R .* Tc);

dTdh = 0.025;
tol = 1e-10;
maxiter = 500;

Cp_coeffs = [
    31.15,     -0.014,    2.68e-5,   -1.17e-8;
    19.79,      0.073,   -5.60e-5,    1.72e-8;
    19.25,      0.052,    1.20e-5,   -1.13e-8;
     5.41,      0.178,   -6.94e-5,    8.71e-9;
    -4.22,      0.306,   -1.59e-4,    3.21e-8;
    -1.39,      0.385,   -1.83e-4,    2.90e-8;
     9.49,      0.331,   -1.11e-4,   -2.82e-9;
    -9.52,      0.507,   -2.73e-4,    5.72e-8;
    -3.63,      0.487,   -2.58e-4,    5.30e-8;
    -4.41,      0.582,   -3.12e-4,    6.49e-8;
    -5.15,      0.676,   -3.65e-4,    7.66e-8;
    -6.10,      0.771,   -4.20e-4,    8.85e-8;
     3.14,      0.677,   -1.93e-4,   -2.98e-8;
    25.20,      0.830,   -3.23e-4,    4.06e-8;
    30.14,      0.993,   -3.87e-4,    4.86e-8;
    36.83,      1.213,   -4.72e-4,    5.93e-8;
    43.82,      1.443,   -5.62e-4,    7.06e-8;
    49.52,      1.630,   -6.35e-4,    7.98e-8;
    57.03,      1.878,   -7.31e-4,    9.19e-8;
    66.63,      2.194,   -8.55e-4,    1.07e-7;
    82.18,      2.706,   -1.05e-3,    1.32e-7;
   115.26,      3.795,   -1.48e-3,    1.86e-7];

% H_ig_ref = [8330.789; 19459.101; 2.642; 9761.134; 19519.622; 29278.121; 29278.121; 39036.609; 39036.609; 48795.1026; ...
%             58553.595; 68312.084; 78069.882; 86301.327; 105418.986; 131284.869; 158312.556; 180345.162; 209390.365; ...
%             246519.548; 306655.264; 434614.185];

 H_ig_ref = [8008.5; 19023.5; 0.0; 10013.1; 20020.1; 28984.2; 28984.2; 39079.8;29054.7; 49022.7; 55626.4; 63278.8; 75018.2; 86305.5; 105436.4;131245.8; 158307.4; 180360.2; 209369.7; 246518.5; 306624.4; 434697.4];

idx_C1  = 3;
idx_C7p = 11:22;

h_oil_base = [175, 204, 228, 327];
h_gas_base = [120, 100, 80, 60, 40, 20, 0];

vt_ids    = [0, 5];
vt_names  = {'PR', 'PR + Abudour VT'};
n_vt = length(vt_ids);

results = struct();

for m = 1:n_vt
    vt_method = vt_ids(m);
    fprintf('\n=== Running Haase + %s (VT=%d) ===\n', vt_names{m}, vt_method);

    vt_params = struct();
    vt_params.Vc = Vc;
    vt_params.components = components;
    vt_params.Zc = Zc;
    if vt_method == 6
        vt_params.c_custom = c_JY;
    end

    [GOC_h, ~, ~, GOC_T] = detect_sgoc(1, comp_ref, press_ref, temp_ref, h_ref, dTdh, [80 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    fprintf('  GOC = %.3f m\n', GOC_h);

    h_oil = [GOC_h, GOC_h+2, h_oil_base];
    n_oil = length(h_oil);
    oil_P = zeros(n_oil,1); oil_Pbub = zeros(n_oil,1);
    oil_rho = zeros(n_oil,1); oil_comp = zeros(n_oil, n);
    oil_GOR = zeros(n_oil,1);
    oil_Bo = zeros(n_oil,1); oil_Bg = zeros(n_oil,1);

    for i = 1:n_oil
        [ch,Ph,Th,Pb,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        oil_P(i) = Ph; oil_Pbub(i) = Pb; oil_comp(i,:) = ch';
        try, oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, oil_rho(i) = NaN; end
        try, [oil_GOR(i), oil_Bo(i), oil_Bg(i)] = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, 'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'PR'); catch ME, fprintf('  GOR failed at h=%.1f: %s\n', h_oil(i), ME.message); oil_GOR(i) = NaN; oil_Bo(i) = NaN; oil_Bg(i) = NaN; end
    end

    [~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter);

    h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
    n_gas = length(h_gas);
    gas_P = zeros(n_gas,1); gas_Pdew = zeros(n_gas,1);
    gas_rho = zeros(n_gas,1); gas_comp = zeros(n_gas, n);
    gas_GOR = zeros(n_gas,1);
    gas_Bo = zeros(n_gas,1); gas_Bg = zeros(n_gas,1);

    for i = 1:n_gas
        [ch,Ph,Th,~,Pd] = main_hasse(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T,0.025, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        gas_P(i) = Ph; gas_Pdew(i) = Pd; gas_comp(i,:) = ch';
        try, gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, gas_rho(i) = NaN; end
        try, [gas_GOR(i), gas_Bo(i), gas_Bg(i)] = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, 'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'PR', 'verbose', true); catch ME, fprintf('  GOR failed at h=%.1f: %s\n', h_gas(i), ME.message); gas_GOR(i) = NaN; gas_Bo(i) = NaN; gas_Bg(i) = NaN; end
    end

    all_h    = [flipud(h_gas'); h_oil'];
    all_comp = [flipud(gas_comp); oil_comp];

    results(m).name = vt_names{m};
    results(m).GOC  = GOC_h;
    results(m).h    = all_h;
    results(m).P    = [flipud(gas_P);  oil_P];
    results(m).Psat = [flipud(gas_Pdew); oil_Pbub];
    results(m).rho  = [flipud(gas_rho); oil_rho];
    results(m).C1   = all_comp(:, idx_C1) * 100;
    results(m).C7p  = sum(all_comp(:, idx_C7p), 2) * 100;
    results(m).GOR  = [flipud(gas_GOR); oil_GOR];
    results(m).Bo   = [flipud(gas_Bo); oil_Bo];
    results(m).Bg   = [flipud(gas_Bg); oil_Bg];
end

%% Experimental data — Pedersen (2015), Reservoir 1
exp_h    = [0; 175; 204; 228; 327];
exp_P    = [279; 284; 286; 287; 293];
exp_Psat = [270; 272; 267; 265; 242];
exp_C1   = [75.66; 50.04; 49.88; 48.89; 45.66];
exp_C7p  = [1.29+1.06+0.54+1.57; 3.50+3.75+2.28+15.88; 3.49+3.73+2.31+16.11; 3.62+3.85+2.36+16.70; 4.01+4.28+2.58+17.66];
exp_GOR  = [3739; 281.5; 280.9; 264.5; 227.8];


GOC_field = 140;

%% ================================================================
%  STYLE — individual thesis figures
%% ================================================================
colors = {[0.00 0.00 0.00], ...     % black   — PR (no VT)
          [0.00 0.45 0.74]};        % blue    — Abudour

lstyles = {'-', '--'};
lw       = 2.0;
lw_sat   = 1.4;
mk_sz    = 7;
fn       = 'Times New Roman';

ax_pos = [0.16 0.16 0.80 0.80];
y_lim  = [0 340];

fmt_ax = @(ax) set(ax, 'YDir', 'reverse', ...
    'FontName', fn, 'FontSize', 11, ...
    'TickDir', 'in', 'TickLength', [0.015 0.015], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 0.8);

fmt_legend = @(lg) set(lg, 'FontSize', 9, 'FontName', fn, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], 'Color', 'w', ...
    'Interpreter', 'tex');

leg_str = cell(1, n_vt);
for m = 1:n_vt
    leg_str{m} = sprintf('%s (GOC %.0f m)', results(m).name, results(m).GOC);
end

%% ===== Figure 1: Pressure & Psat =====
fig1 = figure('Units','centimeters','Position',[2 2 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    h_lines(m) = plot(results(m).P/1e5, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
    v = ~isnan(results(m).Psat) & results(m).Psat > 0;
    plot(results(m).Psat(v)/1e5, results(m).h(v), lstyles{m}, 'Color', colors{m}, 'LineWidth', lw_sat, 'HandleVisibility','off');
end
h_psat_dummy = plot(NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', lw_sat);
h_exp_P    = plot(exp_P,    exp_h, 'ko', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
h_exp_Psat = plot(exp_Psat, exp_h, 'k^', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

fmt_ax(ax); ylim(y_lim); xlim([215 305]);
xlabel('Pressure (bar)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',       'FontSize', 12, 'FontName', fn);

lg = legend([h_exp_P; h_exp_Psat; h_field_GOC; h_psat_dummy; h_lines], ...
    [{'Measured P', 'Measured P_{sat}', sprintf('Field GOC (%d m)', GOC_field), ...
      }, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig1, 'Res1_Fig_pressure', '-dpdf', '-r600');
print(fig1, 'Res1_Fig_pressure', '-dpng', '-r600');

%% ===== Figure 2: Density =====
fig2 = figure('Units','centimeters','Position',[4 3 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    v = ~isnan(results(m).rho) & results(m).rho > 0;
    h_lines(m) = plot(results(m).rho(v), results(m).h(v), lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
exp_rho  = [285.4; 628.6; 650.8; 670.5; 720.7];
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);
h_exp_P    = plot(exp_rho,    exp_h, 'ko', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(ax); ylim(y_lim);
xlabel('Density (kg/m^{3})', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',           'FontSize', 12, 'FontName', fn);

lg = legend([h_field_GOC; h_lines], [{sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig2, 'Res1_Fig_density', '-dpdf', '-r600');
print(fig2, 'Res1_Fig_density', '-dpng', '-r600');

%% ===== Figure 3: C1 mol% =====
fig3 = figure('Units','centimeters','Position',[6 4 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    h_lines(m) = plot(results(m).C1, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
h_exp = plot(exp_C1, exp_h, 'ko', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

fmt_ax(ax); ylim(y_lim);
xlabel('C_{1} (mol%)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',     'FontSize', 12, 'FontName', fn);

lg = legend([h_exp; h_field_GOC; h_lines], [{'Measured', sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig3, 'Res1_Fig_C1', '-dpdf', '-r600');
print(fig3, 'Res1_Fig_C1', '-dpng', '-r600');

%% ===== Figure 4: C7+ mol% =====
fig4 = figure('Units','centimeters','Position',[8 5 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    h_lines(m) = plot(results(m).C7p, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
h_exp = plot(exp_C7p, exp_h, 'ko', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

fmt_ax(ax); ylim(y_lim);
xlabel('C_{7+} (mol%)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',      'FontSize', 12, 'FontName', fn);

lg = legend([h_exp; h_field_GOC; h_lines], [{'Measured', sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig4, 'Res1_Fig_C7p', '-dpdf', '-r600');
print(fig4, 'Res1_Fig_C7p', '-dpng', '-r600');

%% ===== Figure 5: GOR (linear scale) =====
fig5 = figure('Units','centimeters','Position',[10 6 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    v = ~isnan(results(m).GOR) & results(m).GOR > 0;
    h_lines(m) = plot(results(m).GOR(v), results(m).h(v), lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
v_exp = ~isnan(exp_GOR);
if any(v_exp)
    h_exp = plot(exp_GOR(v_exp), exp_h(v_exp), 'ko', 'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
end
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

fmt_ax(ax); ylim(y_lim);
set(ax, 'XScale', 'log');
xlabel('GOR (Sm^{3}/Sm^{3})', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',             'FontSize', 12, 'FontName', fn);

if any(v_exp)
    lg = legend([h_exp; h_field_GOC; h_lines], [{'Measured', sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
else
    lg = legend([h_field_GOC; h_lines], [{sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
end
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig5, 'Res1_Fig_GOR', '-dpdf', '-r600');
print(fig5, 'Res1_Fig_GOR', '-dpng', '-r600');

%% ===== Figure 6: Bo (oil zone only, h > GOC+1) =====
fig6 = figure('Units','centimeters','Position',[12 7 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

all_Bo = [];
h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    v = ~isnan(results(m).Bo) & results(m).Bo > 0 & results(m).h > results(m).GOC + 1;
    if any(v)
        h_lines(m) = plot(results(m).Bo(v), results(m).h(v), [lstyles{m} 'o'], 'Color', colors{m}, ...
            'LineWidth', lw, 'MarkerSize', 4, 'MarkerFaceColor', colors{m});
        all_Bo = [all_Bo; results(m).Bo(v)];
    else
        h_lines(m) = plot(NaN, NaN, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
    end
end
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

Bo_max = max(all_Bo);
y_min_Bo = floor(min(arrayfun(@(s) s.GOC, results)));
y_max_Bo = 340;

fmt_ax(ax); ylim([y_min_Bo y_max_Bo]); xlim([0 Bo_max*1.10]);
xlabel('B_{o} (rm^{3}/Sm^{3})', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',               'FontSize', 12, 'FontName', fn);

lg = legend([h_field_GOC; h_lines], [{sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig6, 'Res1_Fig_Bo', '-dpdf', '-r600');
print(fig6, 'Res1_Fig_Bo', '-dpng', '-r600');

%% ===== Figure 7: Bg (gas zone only, h < GOC-1) =====
fig7 = figure('Units','centimeters','Position',[14 8 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

all_Bg = [];
h_lines = gobjects(n_vt, 1);
for m = 1:n_vt
    v = ~isnan(results(m).Bg) & results(m).Bg > 0 & results(m).h < results(m).GOC - 1;
    if any(v)
        h_lines(m) = plot(results(m).Bg(v)*1e3, results(m).h(v), [lstyles{m} 'o'], 'Color', colors{m}, ...
            'LineWidth', lw, 'MarkerSize', 4, 'MarkerFaceColor', colors{m});
        all_Bg = [all_Bg; results(m).Bg(v)*1e3];
    else
        h_lines(m) = plot(NaN, NaN, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
    end
end
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

Bg_max = max(all_Bg);
y_min_Bg = 0;
y_max_Bg = ceil(max(arrayfun(@(s) s.GOC, results)));

fmt_ax(ax); ylim([y_min_Bg y_max_Bg]); xlim([0 Bg_max*1.10]);
xlabel('B_{g} (\times 10^{-3} rm^{3}/Sm^{3})', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',                              'FontSize', 12, 'FontName', fn);

lg = legend([h_field_GOC; h_lines], [{sprintf('Field GOC (%d m)', GOC_field)}, leg_str], 'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

print(fig7, 'Res1_Fig_Bg', '-dpdf', '-r600');
print(fig7, 'Res1_Fig_Bg', '-dpng', '-r600');

fprintf('\nAll seven figures saved.\n');
fprintf('GOC summary:\n');
for m = 1:n_vt
    fprintf('  %-22s  %.0f m\n', results(m).name, results(m).GOC);
end
fprintf('  Field                   ~140 m\n');