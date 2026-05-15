%% Demo: Isothermal Compositional Grading — VT Comparison
%  12-component fluid, PR EOS, main_hasse with dTdh = 0
%  VT = 0 (no VT), 5 (Abudour)
%%
clear; clc; close all;

components = {'CO2','N2','C1','C2','C3','iC4','nC4','iC5','nC5','nC6','QnC11','QnC33'};
n = 12;

comp_ref = [0.12; 0.13; 54.98; 5.44; 4.88; 1.88; 2.48; 1.45; 1.30; 1.91; 21.53; 3.90];
comp_ref = comp_ref / sum(comp_ref);

Pc = [73.76; 33.94; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; 29.69; 23.43; 9.90]*1e5;
Tc = [304.2; 126.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.4; 642.51; 929.34];
acentric = [0.2250; 0.0400; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.2960; 0.4882; 1.0854];
M_gmol = [44.01; 28.01; 16.04; 30.07; 44.10; 58.12; 58.12; 72.15; 72.15; 86.18; 151.16; 459.34];
Zc = 0.29056 - 0.08775 .* acentric;

BIP = zeros(n);
BIP(3,4:end) = [0.0256, 0.0270, 0.0284, 0.0284, 0.0298, 0.0298, 0.0312, 0.0486, 0.085];
BIP(4:end,3) = BIP(3,4:end);

h_ref = 100;
press_ref = 293.02e5;
temp_ref = 361;
dTdh = 0;
eos_type = 'PR';
tol = 1e-10;
maxiter = 500;

Cp_coeffs = zeros(n, 4);
H_ig_ref = zeros(n, 1);

depth_range = [0, 120];

%% ===== VT method definitions =====
R = 8.3144598;

vt_methods = [0, 5];
vt_labels  = {'PR (no VT)', 'PR + Abudour'};
n_vt = length(vt_methods);

results = struct();

for iv = 1:n_vt
    vt_method = vt_methods(iv);
    fprintf('\n========== VT = %d (%s) ==========\n', vt_method, vt_labels{iv});

    vt_params = struct();
    vt_params.components = components;
    vt_params.Zc = Zc;
    vt_params.Vc = Zc .* R .* Tc ./ Pc * 1e6;

    [GOC_depth, GOC_pressure, ~, GOC_temp] = detect_sgoc(1, comp_ref, press_ref, temp_ref, ...
        h_ref, dTdh, depth_range, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    fprintf('GOC = %.2f m (internal), %.2f m (above datum)\n', GOC_depth, h_ref - GOC_depth);

    h_oil = unique([GOC_depth, ceil(GOC_depth):1:depth_range(2)], 'sorted')';
    n_oil = length(h_oil);
    comp_oil = zeros(n_oil, n); P_oil = zeros(n_oil, 1);
    Pbub_oil = NaN(n_oil, 1); rho_oil = NaN(n_oil, 1);

    for i = 1:n_oil
        try
            [z_h, P_h, T_h, ~, ~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, ...
                dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
            comp_oil(i,:) = z_h(:)'; P_oil(i) = P_h;
            Pb = calc_pbub(z_h, P_h, T_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
            if ~isnan(Pb) && Pb <= P_h, Pbub_oil(i) = Pb; end
            try rho_oil(i) = calculate_density(z_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, end
        catch, end
    end

    [~, y_eq] = pressbub_multicomp_newton(comp_oil(1,:)', P_oil(1), temp_ref, ...
        Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
    y_eq = max(y_eq(:), 1e-15); y_eq = y_eq / sum(y_eq);

    h_gas = unique([depth_range(1):1:floor(GOC_depth), GOC_depth], 'sorted')';
    n_gas = length(h_gas);
    comp_gas = zeros(n_gas, n); P_gas = zeros(n_gas, 1);
    Pdew_gas = NaN(n_gas, 1); rho_gas = NaN(n_gas, 1);

    for i = 1:n_gas
        try
            [z_h, P_h, T_h, ~, ~] = main_hasse(h_gas(i), GOC_depth, y_eq, GOC_pressure, temp_ref, ...
                dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
            comp_gas(i,:) = z_h(:)'; P_gas(i) = P_h;
            Pd = calc_pdew(z_h, P_h, T_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
            if ~isnan(Pd) && Pd <= P_h, Pdew_gas(i) = Pd; end
            try rho_gas(i) = calculate_density(z_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, end
        catch, end
    end

    results(iv).vt       = vt_method;
    results(iv).label    = vt_labels{iv};
    results(iv).GOC      = h_ref - GOC_depth;
    results(iv).h_plot   = h_ref - [h_gas; h_oil];
    results(iv).P        = [P_gas; P_oil];
    results(iv).Psat     = [Pdew_gas; Pbub_oil];
    results(iv).C1       = [comp_gas(:,3); comp_oil(:,3)];
    results(iv).C7p      = [sum(comp_gas(:,11:12),2); sum(comp_oil(:,11:12),2)];
    results(iv).rho      = [rho_gas; rho_oil];
end

%% ---- Schulte (1980) SPE 9235, Table 3 data ----
schulte_oil_h_plot = [0; 22];
schulte_oil_P = [293.02; 291.84] * 1e5;
schulte_oil_Pbub = [290.14; 291.84] * 1e5;
schulte_oil_C1 = [0.5498; 0.5559];
schulte_oil_C7p = [0.0294+0.1216+0.1033; 0.0288+0.1203+0.0987];

schulte_gas_h_plot = [22; 50; 100];
schulte_gas_P = [291.84; 291.23; 290.13] * 1e5;
schulte_gas_Pdew = [291.84; 289.16; 284.61] * 1e5;
schulte_gas_C1 = [0.8415; 0.8429; 0.8454];
schulte_gas_C7p = [0.0029+0.0087+0.0058; 0.0028+0.0084+0.0057; 0.0027+0.0077+0.0054];

schulte_h_plot = [schulte_gas_h_plot; schulte_oil_h_plot];
schulte_P      = [schulte_gas_P; schulte_oil_P];
schulte_Psat   = [schulte_gas_Pdew; schulte_oil_Pbub];
schulte_C1     = [schulte_gas_C1; schulte_oil_C1];
schulte_C7p    = [schulte_gas_C7p; schulte_oil_C7p];

%% ===== Style =====
col_vt = [
    0.00 0.45 0.74      % blue   — no VT
    0.85 0.33 0.10      % orange — Abudour
];
col_sch = [0.00 0.00 0.00];

lw     = 2.0;
lw_sat = 1.4;
mk_sz  = 7;

ax_pos = [0.15 0.16 0.81 0.80];

%% ===== Figure 1: Pressure & Psat =====
fig1 = figure('Units','centimeters','Position',[2 2 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_vt, 1);
for iv = 1:n_vt
    h_lines(iv) = plot(results(iv).P/1e5, results(iv).h_plot, '-', ...
        'Color', col_vt(iv,:), 'LineWidth', lw);
    v = ~isnan(results(iv).Psat) & results(iv).Psat > 0;
    plot(results(iv).Psat(v)/1e5, results(iv).h_plot(v), '--', ...
        'Color', col_vt(iv,:), 'LineWidth', lw_sat, 'HandleVisibility', 'off');
end
h_psat_dummy = plot(NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', lw_sat);
h_exp_P    = plot(schulte_P/1e5, schulte_h_plot, 's', 'Color', col_sch, ...
    'MarkerFaceColor', col_sch, 'MarkerSize', mk_sz);
h_exp_Psat = plot(schulte_Psat/1e5, schulte_h_plot, 'd', 'Color', col_sch, ...
    'MarkerSize', mk_sz, 'LineWidth', 1.0);

set(ax, 'FontName','Times New Roman', 'FontSize', 11, ...
        'TickDir','in', 'TickLength', [0.015 0.015], ...
        'XMinorTick','on', 'YMinorTick','on', 'LineWidth', 0.8);
xlabel('Pressure (bar)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Height above datum (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

lg = legend([h_lines; h_psat_dummy; h_exp_P; h_exp_Psat], ...
    [vt_labels, {'P_{sat} (model)', 'Schulte P', 'Schulte P_{sat}'}], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ===== Figure 2: Composition (C1 & C7+) =====
fig2 = figure('Units','centimeters','Position',[4 3 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_C1 = gobjects(n_vt, 1);
for iv = 1:n_vt
    h_C1(iv) = plot(results(iv).C1, results(iv).h_plot, '-', ...
        'Color', col_vt(iv,:), 'LineWidth', lw);
    plot(results(iv).C7p, results(iv).h_plot, '--', ...
        'Color', col_vt(iv,:), 'LineWidth', lw, 'HandleVisibility', 'off');
end
h_C7p_dummy = plot(NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', lw);
h_exp_C1   = plot(schulte_C1, schulte_h_plot, 's', 'Color', col_sch, ...
    'MarkerFaceColor', col_sch, 'MarkerSize', mk_sz);
h_exp_C7p  = plot(schulte_C7p, schulte_h_plot, 'd', 'Color', col_sch, ...
    'MarkerSize', mk_sz, 'LineWidth', 1.0);

set(ax, 'FontName','Times New Roman', 'FontSize', 11, ...
        'TickDir','in', 'TickLength', [0.015 0.015], ...
        'XMinorTick','on', 'YMinorTick','on', 'LineWidth', 0.8);
xlabel('Mole fraction', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Height above datum (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

lg = legend([h_C1; h_C7p_dummy; h_exp_C1; h_exp_C7p], ...
    [strcat(vt_labels, {' (C_{1})'}), {'C_{7+} (dashed)', 'Schulte C_{1}', 'Schulte C_{7+}'}], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ===== Figure 3: Density =====
fig3 = figure('Units','centimeters','Position',[6 4 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_rho = gobjects(n_vt, 1);
for iv = 1:n_vt
    v = ~isnan(results(iv).rho) & results(iv).rho > 0;
    h_rho(iv) = plot(results(iv).rho(v), results(iv).h_plot(v), '-', ...
        'Color', col_vt(iv,:), 'LineWidth', lw);
end

set(ax, 'FontName','Times New Roman', 'FontSize', 11, ...
        'TickDir','in', 'TickLength', [0.015 0.015], ...
        'XMinorTick','on', 'YMinorTick','on', 'LineWidth', 0.8);
xlabel('Density (kg/m^{3})', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Height above datum (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

lg = legend(h_rho, vt_labels, ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ========== LOCAL HELPERS ==========

function Pb = calc_pbub(z, P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp)
    Pb = NaN;
    try
        [Pb, ~] = pressbub_multicomp_newton(z(:), P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp);
        if ~isreal(Pb) || Pb <= 0 || ~isfinite(Pb), Pb = NaN; end
    catch
    end
    if isnan(Pb)
        try
            [Pb, ~] = pressbub_multicomp_ss(z(:), P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp);
            if ~isreal(Pb) || Pb <= 0 || ~isfinite(Pb), Pb = NaN; end
        catch
            Pb = NaN;
        end
    end
end

function Pd = calc_pdew(z, P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp)
    Pd = NaN;
    try
        [Pd, ~] = pressdew_multicomp_newton(z(:), P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp);
        if ~isreal(Pd) || Pd <= 0 || ~isfinite(Pd), Pd = NaN; end
    catch
    end
    if isnan(Pd)
        try
            [Pd, ~] = pressdew_multicomp_ss(z(:), P, T, Pc, Tc, w, BIP, tol, maxiter, eos, comp);
            if ~isreal(Pd) || Pd <= 0 || ~isfinite(Pd), Pd = NaN; end
        catch
            Pd = NaN;
        end
    end
end