clear; clc; close all;

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','nC6','nC7','nC8','nC9','nC10','nC14','nC17','nC20','nC23','nC28','nC34','nC39-48','C49-C80'};
comp_ref = [0.395; 2.060; 53.8715362; 7.589; 5.575; 1.009; 2.514; 0.900; 1.396; 1.557; ...
            2.630; 2.823; 1.897; 4.406; 2.479; 1.941; 1.520; 1.526; 1.325; 1.145; 0.829; 0.610];
comp_ref = comp_ref ./ sum(comp_ref);

M_gmol = [28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; ...
          96.000; 107.000; 121.000; 152.857; 205.131; 249.627; 289.518; 336.857; 399.723; 481.458; 595.685; 811.198];
Tc = [-146.950; 31.050; -82.550; 32.250; 96.650; 134.950; 152.050; 187.250; 196.450; 234.250; ...
      262.184; 282.716; 304.023; 346.593; 401.967; 442.808; 476.719; 514.042; 560.087; 616.024; 689.007; 821.700] + 273.15;
Pc = [33.958; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; 29.69; ...
      31.95; 29.76; 26.67; 22.04; 18.06; 16.32; 15.43; 14.71; 14.12; 13.69; 13.43; 13.40] * 1e5;
acentric = [0.0377; 0.2250; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.2960; ...
            0.4679; 0.4999; 0.5399; 0.6321; 0.7673; 0.8746; 0.9642; 1.0608; 1.1710; 1.2799; 1.3595; 1.2257];

n = 22;
h_ref = 3682.8;
P_ref = 380.20e5;
T_ref = 138.70 + 273.15;
dTdh  = 0.0260;
R = 8.3144598;

Zc = 0.29056 - 0.08775 .* acentric;
for i = 2:length(Zc)
    if Zc(i) >= Zc(i-1)
        Zc(i) = Zc(i-1) - 0.008;
    end
end
Vc = Zc .* R .* Tc ./ Pc;

c_custom = [0.92; 3.03; 0.63; 2.63; 5.06; 7.29; 7.86; 10.93; 12.18; 17.98; ...
            9.65; 13.90; 20.99; 34.32; 45.88; 45.87; 39.23; 26.41; 3.39; -33.93; -94.44; -219.18];

c_peneloux = [-4.23; -1.64; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12; 1.39; ...
               6.70; 10.80; 14.52; 21.55; 26.56; 25.75; 22.30; 15.13; 1.78; -19.81; -55.16; -124.90];

BIP = zeros(n);
BIP(2,1) = -0.0315; BIP(3,1) = 0.0278; BIP(4,1) = 0.0407; BIP(5,1) = 0.0763;
BIP(6,1) = 0.0944; BIP(7,1) = 0.0700; BIP(8,1) = 0.0867; BIP(9,1) = 0.0878; BIP(10,1) = 0.0800;
for i = 11:22, BIP(i,1) = 0.0800; end
for i = 3:10,  BIP(i,2) = 0.1200; end
for i = 11:22, BIP(i,2) = 0.1000; end
BIP = BIP + BIP';

Cp_coeffs = [
    31.15,   -0.014,    2.68e-5,  -1.17e-8;
    19.79,    0.073,   -5.60e-5,   1.72e-8;
    19.25,    0.052,    1.20e-5,  -1.13e-8;
     5.41,    0.178,   -6.94e-5,   8.71e-9;
    -4.22,    0.306,   -1.59e-4,   3.21e-8;
    -1.39,    0.385,   -1.85e-4,   2.90e-8;
     9.49,    0.331,   -1.11e-4,  -2.82e-9;
     9.52,    0.507,   -2.73e-4,   5.72e-8;
    -3.63,    0.487,   -2.58e-4,   5.30e-8;
    -4.41,    0.582,   -3.12e-4,   6.49e-8;
     9.58,    0.576,   -2.05e-4,   0.000;
    -4.91,    0.612,   -2.33e-4,   0.000;
    -1.87,    0.680,   -2.69e-4,   0.000;
     2.31,    0.859,   -3.48e-4,   0.000;
     5.06,    1.151,   -4.67e-4,   0.000;
     7.19,    1.402,   -5.66e-4,   0.000;
     8.41,    1.629,   -6.57e-4,   0.000;
    10.06,    1.899,   -7.65e-4,   0.000;
    12.3129,  2.258,   -9.08e-4,   0.000;
    14.70,    2.721,   -1.09e-3,   0.000;
    18.27,    3.369,   -1.36e-3,   0.000;
    26.84,    4.652,   -1.87e-3,   0.000];

H_ig_ref = [8330.8; 19459.1; 2.6; 9761.1; 19519.6; 29278.1; 29278.1; 39036.6; 39036.6; 48795.1; ...
            55628.2; 63280.8; 73020.5; 95183.2; 131550.0; 162505.2; 190257.4; 223190.8; 266926.5; 323789.1; 403255.9; 553186.7];

cases = struct();

cases(1).name = 'SRK + no VT';
cases(1).vt   = 7;
cases(1).params = struct('c_custom', c_peneloux, 'Vc', Vc, 'components', {components}, 'Zc', Zc);
cases(1).hardcode_GOC = NaN;

cases(2).name = 'SRK + Chen-Li';
cases(2).vt   = 9;
cases(2).params = struct('c_custom', c_custom, 'Vc', Vc, 'components', {components}, 'Zc', Zc);
cases(2).hardcode_GOC = 3646;

cases(3).name = 'SRK + Constant VT';
cases(3).vt   = 12;
cases(3).params = struct('c_custom', c_custom, 'Vc', Vc, 'components', {components}, 'Zc', Zc);
cases(3).hardcode_GOC = NaN;

n_cases = 3;

h_extra_oil = unique([3648; 3651.1; 3655; 3661.6; 3670; 3676.0; ...
                      linspace(3647, 3680, 40)']);
h_extra_gas = flipud(unique([3644.3; 3642; 3638.2; 3630; ...
                             linspace(3624, 3645.5, 25)']));

comp_gas_GOC = [0.566; 2.279; 66.572; 8.069; 5.434; 0.940; 2.248; 0.766; 1.163; 1.226; ...
                1.809; 1.849; 1.192; 2.518; 1.206; 0.792; 0.512; 0.402; 0.247; 0.135; 0.054; 0.021];
comp_gas_GOC = comp_gas_GOC ./ sum(comp_gas_GOC);

for k = 1:n_cases
    vm = cases(k).vt;
    vp = cases(k).params;
    fprintf('\n========== Case %d: %s (vt=%d) ==========\n', k, cases(k).name, vm);

    if ~isnan(cases(k).hardcode_GOC)
        GOC_h = cases(k).hardcode_GOC;
    else
        [GOC_h, ~, ~, ~] = detect_sgoc(1, comp_ref, P_ref, T_ref, h_ref, dTdh, ...
            [3600 3680], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
    end
    fprintf('  GOC: %.1f m\n', GOC_h);
    cases(k).GOC = GOC_h;
    GOC_T = T_ref + dTdh * (GOC_h - h_ref);

    h_oil     = [GOC_h; h_extra_oil];
    h_gas_arr = [GOC_h; h_extra_gas];
    n_oil = length(h_oil);
    n_gas = length(h_gas_arr);

    vm_gas = vm;
    vp_gas = vp;
    if vm == 9
        vm_gas = 12;
        vp_gas = cases(3).params;
    end

    oil_P = NaN(n_oil,1); oil_T = NaN(n_oil,1);
    oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
    for i = 1:n_oil
        try
            [ch, Ph, Th, ~, ~] = main_hasse(h_oil(i), h_ref, comp_ref, P_ref, T_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
            oil_P(i) = Ph; oil_T(i) = Th; oil_comp(i,:) = ch(:)';
            [rho,~,~,~] = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
            oil_rho(i) = rho;
        catch ME
            fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message);
        end
    end

    P_GOC = oil_P(1);

    gas_P = NaN(n_gas,1); gas_T = NaN(n_gas,1);
    gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
    for i = 1:n_gas
        try
            [ch, Ph, Th, ~, ~] = main_hasse(h_gas_arr(i), GOC_h, comp_gas_GOC, P_GOC, GOC_T, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm_gas, vp_gas);
            gas_P(i) = Ph; gas_T(i) = Th; gas_comp(i,:) = ch(:)';
            [rho,~,~,~] = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm_gas, vp_gas);
            gas_rho(i) = rho;
        catch ME
            fprintf('  Gas error at h=%.1f: %s\n', h_gas_arr(i), ME.message);
        end
    end

    oil_GOR = NaN(n_oil,1);
    gas_GOR = NaN(n_gas,1);
    for i = 1:n_oil
        if all(oil_comp(i,:) > 0) && ~isnan(oil_P(i))
            try
                [GOR_val, ~, ~, ~] = calculate_GOR_STO(oil_comp(i,:)', oil_P(i), oil_T(i), ...
                    Pc, Tc, acentric, BIP, M_gmol, ...
                    'vt_method', vm, 'vt_params', vp, 'eos_type', 'SRK');
                oil_GOR(i) = GOR_val;
            catch ME
                fprintf('  GOR error (oil) at h=%.1f: %s\n', h_oil(i), ME.message);
            end
        end
    end
    for i = 1:n_gas
        if all(gas_comp(i,:) > 0) && ~isnan(gas_P(i))
            try
                [GOR_val, ~, ~, ~] = calculate_GOR_STO(gas_comp(i,:)', gas_P(i), gas_T(i), ...
                    Pc, Tc, acentric, BIP, M_gmol, ...
                    'vt_method', vm_gas, 'vt_params', vp_gas, 'eos_type', 'SRK');
                gas_GOR(i) = GOR_val;
            catch ME
                fprintf('  GOR error (gas) at h=%.1f: %s\n', h_gas_arr(i), ME.message);
            end
        end
    end

    oil_C1  = oil_comp(:,3)  * 100;
    gas_C1  = gas_comp(:,3)  * 100;
    oil_C7p = sum(oil_comp(:,11:end), 2) * 100;
    gas_C7p = sum(gas_comp(:,11:end), 2) * 100;

    cases(k).h_oil = h_oil;         cases(k).h_gas  = h_gas_arr;
    cases(k).oil_rho = oil_rho;     cases(k).gas_rho = gas_rho;
    cases(k).oil_comp = oil_comp;   cases(k).gas_comp = gas_comp;
    cases(k).oil_P = oil_P;         cases(k).gas_P = gas_P;
    cases(k).oil_T = oil_T;         cases(k).gas_T = gas_T;
    cases(k).oil_GOR = oil_GOR;     cases(k).gas_GOR = gas_GOR;
    cases(k).oil_C1  = oil_C1;      cases(k).gas_C1  = gas_C1;
    cases(k).oil_C7p = oil_C7p;     cases(k).gas_C7p = gas_C7p;

    cases(k).all_h    = [flipud(h_gas_arr); h_oil];
    cases(k).all_rho  = [flipud(gas_rho);   oil_rho];
    cases(k).all_GOR  = [flipud(gas_GOR);   oil_GOR];
    cases(k).all_C1   = [flipud(gas_C1);    oil_C1];
    cases(k).all_C7p  = [flipud(gas_C7p);   oil_C7p];
end

%% Experimental data
exp_h_rho = [3638.2; 3638.2; 3644.3; 3644.3; 3651.1; 3651.1; 3661.6; 3661.6; 3676.0; 3676.0; 3682.8; 3682.8];
exp_rho   = [0.367; 0.374; 0.376; 0.380; 0.570; 0.574; 0.581; 0.584; 0.590; 0.602; 0.597; 0.605] * 1000;

exp_h_gor = [3638.2; 3638.2; 3638.2; 3644.3; 3644.3; 3651.1; 3651.1; 3661.6; 3676.0; 3682.8; 3682.8];
exp_gor   = [1010; 1095; 1150; 1060; 1150; 305; 340; 312; 285; 258; 280];

exp_h_c1  = [3638.2; 3638.2; 3644.3; 3644.3; 3651.1; 3651.1; 3661.6; 3661.6; 3676.0; 3676.0; 3682.8; 3682.8];
exp_c1    = [68.0; 69.7; 68.0; 68.54; 56.142; 56.142; 55.261; 55.261; 54.253; 54.253; 53.871; 53.871];

exp_h_c7p = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];
exp_c7p   = [1.482+1.595+1.031+5.231; ...
             1.519+1.610+1.048+5.080; ...
             2.474+2.583+1.695+13.78; ...
             2.520+2.667+1.779+14.491; ...
             2.579+2.777+1.869+15.282; ...
             2.630+2.823+1.897+15.783];

exp_h_avg   = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];
exp_rho_avg = [367;    376;    574;    581;    595;    601];
exp_GOR_avg = [1086;   1105;   323.0;  311.9;  285.2;  268.3];
exp_C1_avg  = [68.861; 68.546; 56.142; 55.261; 54.253; 53.871];

%% ===== Style =====
colors  = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.00 0.60 0.30]};
lstyles = {'-', '--', '-.'};
mk_face = [0.0 0.0 0.0];
fn      = 'Times New Roman';
lw      = 2.0;
mk_sz   = 7;

ax_pos = [0.16 0.16 0.80 0.80];
y_lim  = [3620 3700];

leg_labels_full = [arrayfun(@(k) cases(k).name, 1:n_cases, 'UniformOutput', false), {'Exp. data (Pedersen 2006)'}];

fmt_ax = @(ax) set(ax, 'YDir','reverse', ...
    'FontName', fn, 'FontSize', 11, ...
    'TickDir','in', 'TickLength', [0.015 0.015], ...
    'XMinorTick','on', 'YMinorTick','on', 'LineWidth', 0.8);

fmt_legend = @(lg) set(lg, 'FontSize', 9, 'FontName', fn, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], 'Color', 'w', ...
    'Interpreter', 'tex');

%% ===== Figure 1: Density =====
fig1 = figure('Units','centimeters','Position',[2 2 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_legs = gobjects(n_cases,1);
for k = 1:n_cases
    col = colors{k}; ls = lstyles{k};
    v_oil = ~isnan(cases(k).oil_rho);
    v_gas = ~isnan(cases(k).gas_rho);

    rho_oil = cases(k).oil_rho(v_oil);
    h_oil_v = cases(k).h_oil(v_oil);
    [h_oil_v, sidx] = sort(h_oil_v);
    rho_oil = rho_oil(sidx);
    if length(h_oil_v) >= 2 && h_oil_v(end) < h_ref
        n_use = min(3, length(h_oil_v));
        p_fit = polyfit(h_oil_v(end-n_use+1:end), rho_oil(end-n_use+1:end), 1);
        rho_oil(end+1) = polyval(p_fit, h_ref);
        h_oil_v(end+1) = h_ref;
    end
    h_legs(k) = plot(rho_oil, h_oil_v, ls, 'Color', col, 'LineWidth', lw);

    plot(cases(k).gas_rho(v_gas), cases(k).h_gas(v_gas), ls, ...
        'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    rho_g = cases(k).gas_rho(1); rho_o = cases(k).oil_rho(1);
    if ~isnan(rho_g) && ~isnan(rho_o)
        plot([rho_g, rho_o], [cases(k).GOC, cases(k).GOC], ls, ...
            'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end
h_exp = plot(exp_rho, exp_h_rho, 'ko', 'MarkerSize', mk_sz, ...
    'MarkerFaceColor', mk_face, 'LineWidth', 0.8);
fmt_ax(ax); xlim([300 650]); ylim(y_lim);
xlabel('Density (kg/m^3)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',         'FontSize', 12, 'FontName', fn);
lg = legend([h_exp; h_legs], [leg_labels_full(end), leg_labels_full(1:n_cases)], ...
    'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

%% ===== Figure 2: GOR =====
fig2 = figure('Units','centimeters','Position',[4 3 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_legs = gobjects(n_cases,1);
for k = 1:n_cases
    col = colors{k}; ls = lstyles{k};
    v_oil = ~isnan(cases(k).oil_GOR);
    v_gas = ~isnan(cases(k).gas_GOR);

    gor_oil = cases(k).oil_GOR(v_oil);
    h_oil_v = cases(k).h_oil(v_oil);
    [h_oil_v, sidx] = sort(h_oil_v);
    gor_oil = gor_oil(sidx);
    if length(h_oil_v) >= 2 && h_oil_v(end) < h_ref
        n_use = min(3, length(h_oil_v));
        p_fit = polyfit(h_oil_v(end-n_use+1:end), gor_oil(end-n_use+1:end), 1);
        gor_oil(end+1) = polyval(p_fit, h_ref);
        h_oil_v(end+1) = h_ref;
    end
    h_legs(k) = plot(gor_oil, h_oil_v, ls, 'Color', col, 'LineWidth', lw);

    plot(cases(k).gas_GOR(v_gas), cases(k).h_gas(v_gas), ls, ...
        'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    gor_g = cases(k).gas_GOR(1); gor_o = cases(k).oil_GOR(1);
    if ~isnan(gor_g) && ~isnan(gor_o)
        plot([gor_g, gor_o], [cases(k).GOC, cases(k).GOC], ls, ...
            'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end
h_exp = plot(exp_gor, exp_h_gor, 'ko', 'MarkerSize', mk_sz, ...
    'MarkerFaceColor', mk_face, 'LineWidth', 0.8);
fmt_ax(ax); xlim([200 1200]); ylim(y_lim);
xlabel('GOR (Sm^3/Sm^3)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',         'FontSize', 12, 'FontName', fn);
lg = legend([h_exp; h_legs], [leg_labels_full(end), leg_labels_full(1:n_cases)], ...
    'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

%% ===== Figure 3: C1 mol% =====
fig3 = figure('Units','centimeters','Position',[6 4 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_legs = gobjects(n_cases,1);
for k = 1:n_cases
    col = colors{k}; ls = lstyles{k};
    v_oil = cases(k).oil_C1 > 0;
    v_gas = cases(k).gas_C1 > 0;
    h_legs(k) = plot(cases(k).oil_C1(v_oil), cases(k).h_oil(v_oil), ls, ...
        'Color', col, 'LineWidth', lw);
    plot(cases(k).gas_C1(v_gas), cases(k).h_gas(v_gas), ls, ...
        'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    c1_g = cases(k).gas_C1(1); c1_o = cases(k).oil_C1(1);
    if c1_g > 0 && c1_o > 0
        plot([c1_g, c1_o], [cases(k).GOC, cases(k).GOC], ls, ...
            'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end
h_exp = plot(exp_c1, exp_h_c1, 'ko', 'MarkerSize', mk_sz, ...
    'MarkerFaceColor', mk_face, 'LineWidth', 0.8);
fmt_ax(ax); xlim([50 72]); ylim(y_lim);
xlabel('C_{1} content (mol%)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',                'FontSize', 12, 'FontName', fn);
lg = legend([h_exp; h_legs], [leg_labels_full(end), leg_labels_full(1:n_cases)], ...
    'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

%% ===== Figure 4: C7+ mol% =====
fig4 = figure('Units','centimeters','Position',[8 5 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_legs = gobjects(n_cases,1);
for k = 1:n_cases
    col = colors{k}; ls = lstyles{k};
    v_oil = cases(k).oil_C7p > 0;
    v_gas = cases(k).gas_C7p > 0;
    h_legs(k) = plot(cases(k).oil_C7p(v_oil), cases(k).h_oil(v_oil), ls, ...
        'Color', col, 'LineWidth', lw);
    plot(cases(k).gas_C7p(v_gas), cases(k).h_gas(v_gas), ls, ...
        'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    c7_g = cases(k).gas_C7p(1); c7_o = cases(k).oil_C7p(1);
    if c7_g > 0 && c7_o > 0
        plot([c7_g, c7_o], [cases(k).GOC, cases(k).GOC], ls, ...
            'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end
h_exp_c7p = plot(exp_c7p, exp_h_c7p, 'ko', 'MarkerSize', mk_sz, ...
    'MarkerFaceColor', mk_face, 'LineWidth', 0.8);
fmt_ax(ax); ylim(y_lim);
xlabel('C_{7+} content (mol%)', 'FontSize', 12, 'FontName', fn);
ylabel('Depth (m)',               'FontSize', 12, 'FontName', fn);
lg = legend([h_exp_c7p; h_legs], [leg_labels_full(end), leg_labels_full(1:n_cases)], ...
    'Location', 'best');
fmt_legend(lg); lg.ItemTokenSize = [22, 10];

%% Summary tables
for k = 1:n_cases
    fprintf('\n=== %s (vt=%d, GOC=%.1f m) ===\n', cases(k).name, cases(k).vt, cases(k).GOC);

    fprintf('\n  Density:\n');
    fprintf('  %-10s %-12s %-12s %-10s\n', 'Depth(m)', 'Rho_model', 'Rho_meas', 'Error(%)');
    for j = 1:length(exp_h_avg)
        [~, idx] = min(abs(cases(k).all_h - exp_h_avg(j)));
        if ~isnan(cases(k).all_rho(idx))
            err = (cases(k).all_rho(idx) - exp_rho_avg(j)) / exp_rho_avg(j) * 100;
            fprintf('  %-10.1f %-12.1f %-12.1f %-10.1f\n', exp_h_avg(j), cases(k).all_rho(idx), exp_rho_avg(j), err);
        end
    end

    fprintf('\n  GOR:\n');
    fprintf('  %-10s %-12s %-12s %-10s\n', 'Depth(m)', 'GOR_model', 'GOR_meas', 'Error(%)');
    for j = 1:length(exp_h_avg)
        [~, idx] = min(abs(cases(k).all_h - exp_h_avg(j)));
        if ~isnan(cases(k).all_GOR(idx))
            err = (cases(k).all_GOR(idx) - exp_GOR_avg(j)) / exp_GOR_avg(j) * 100;
            fprintf('  %-10.1f %-12.1f %-12.1f %-10.1f\n', exp_h_avg(j), cases(k).all_GOR(idx), exp_GOR_avg(j), err);
        end
    end

    fprintf('\n  C1 mol%%:\n');
    fprintf('  %-10s %-12s %-12s %-10s\n', 'Depth(m)', 'C1_model', 'C1_meas', 'Error(%)');
    for j = 1:length(exp_h_avg)
        [~, idx] = min(abs(cases(k).all_h - exp_h_avg(j)));
        if cases(k).all_C1(idx) > 0
            err = (cases(k).all_C1(idx) - exp_C1_avg(j)) / exp_C1_avg(j) * 100;
            fprintf('  %-10.1f %-12.2f %-12.2f %-10.1f\n', exp_h_avg(j), cases(k).all_C1(idx), exp_C1_avg(j), err);
        end
    end
end