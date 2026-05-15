clear; clc; close all;

%% ========================================================================
%  PEDERSEN 2015 RESERVOIR 1 (PR-EOS, No VT) - 5 MODEL COMPARISON
%  ========================================================================
components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7','C8','C9','C10-C11','C12-C13','C14-C16','C17-C18','C19-C21','C22-C24','C25-C29','C30-C37','C38-C80'};

comp_ref = [0.42; 0.71; 49.885; 7.771; 6.741; 1.030; 3.21; 1.16; 1.55; 1.88; ...
            3.490; 3.730; 2.31; 3.293; 2.620; 2.961; 1.480; 1.673; 1.187; 1.264; 0.982; 0.652];
comp_ref = comp_ref ./ sum(comp_ref);

M_gmol = [28.01; 44.01; 16.04; 30.07; 44.10; 58.12; 58.12; 72.15; 72.15; 86.18; ...
          96; 107; 121; 140.09; 167.57; 204.75; 243.6; 275.27; 317.02; 370.39; 456.83; 640.76];
acentric = [0.0400; 0.2250; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.2960; ...
            0.3375; 0.3743; 0.4204; 0.4833; 0.5703; 0.6847; 0.7947; 0.8805; 0.9836; 1.0983; 1.2292; 1.1596];
Tc = [-146.95; 31.05; -82.55; 32.25; 96.65; 134.95; 152.05; 187.25; 196.45; 234.25; ...
      301.63; 323.97; 349.62; 381.84; 422.96; 473.01; 519.43; 555.64; 600.26; 654.86; 737.85; 912.38] + 273.15;
Pc = [33.94; 73.76; 46.00; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; ...
      29.59; 27.37; 25.01; 22.59; 20.12; 17.91; 16.4; 15.55; 14.71; 13.96; 13.16; 12.27] * 1e5;
c_JY = [-4.23; -1.64; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12; 1.39; ...
        15.53; 19.96; 25.14; 31.89; 39.65; 46.52; 49.62; 50.37; 48.37; 42.77; 27.92; -6.04];

R = 8.3144598;
n = length(comp_ref);
h_ref = 204;
press_ref = 286e5;
temp_ref = 94 + 273.15;
dTdh = 0.025;
tau = 4.0;
tol = 1e-10;
maxiter = 500;

BIP = zeros(n);
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

Vc = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0; ...
      432; 492; 548; 600.22; 697.08; 810.90; 910.91; 980.42; 1061.51; 1152.43; 1287.26; 1648.29];
Zc = Pc .* (Vc * 1e-6) ./ (R .* Tc);

Cp_coeffs = [
    31.15,   -0.014,    2.68e-5,  -1.17e-8;
    19.79,    0.073,   -5.60e-5,   1.72e-8;
    19.25,    0.052,    1.20e-5,  -1.13e-8;
     5.41,    0.178,   -6.94e-5,   8.71e-9;
    -4.22,    0.306,   -1.59e-4,   3.21e-8;
    -1.39,    0.385,   -1.83e-4,   2.90e-8;
     9.49,    0.331,   -1.11e-4,  -2.82e-9;
    -9.52,    0.507,   -2.73e-4,   5.72e-8;
    -3.63,    0.487,   -2.58e-4,   5.30e-8;
    -4.41,    0.582,   -3.12e-4,   6.49e-8;
    -5.15,    0.676,   -3.65e-4,   7.66e-8;
    -6.10,    0.771,   -4.20e-4,   8.85e-8;
     3.14,    0.677,   -1.93e-4,  -2.98e-8;
    25.20,    0.830,   -3.23e-4,   4.06e-8;
    30.14,    0.993,   -3.87e-4,   4.86e-8;
    36.83,    1.213,   -4.72e-4,   5.93e-8;
    43.82,    1.443,   -5.62e-4,   7.06e-8;
    49.52,    1.630,   -6.35e-4,   7.98e-8;
    57.03,    1.878,   -7.31e-4,   9.19e-8;
    66.63,    2.194,   -8.55e-4,   1.07e-7;
    82.18,    2.706,   -1.05e-3,   1.32e-7;
   115.26,    3.795,   -1.48e-3,   1.86e-7];

H_ig_ref = [8330.789; 19459.101; 2.642; 9761.134; 19519.622; 29278.121; 29278.121; 39036.609; 39036.609; 48795.1026; ...
            58553.595; 68312.084; 78069.882; 86301.327; 105418.986; 131284.869; 158312.556; 180345.162; 209390.365; ...
            246519.548; 306655.264; 434614.185];

vt_params.c_custom = c_JY;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;
vt_method = 0;

idx_C1  = 3;
idx_C7p = 11:22;

% SR param-free tau
Tc_mix = sum(comp_ref .* Tc);
Tr_mix = temp_ref / Tc_mix;
tau_sr = Tr_mix^0.25 + log(M_gmol) .* acentric.^2 - acentric;
tau_sr = max(min(tau_sr, 6.0), 2.5);

fprintf('=== SR param-free tau ===\n');
fprintf('  Tr_mix = %.4f,  Tr^(1/4) = %.4f\n', Tr_mix, Tr_mix^0.25);
fprintf('  %-10s  %8s  %6s  %8s  %8s\n', 'Component', 'M', 'omega', 'tau=4', 'tau_SR');
fprintf('  %s\n', repmat('-', 1, 48));
for ii = 1:n
    fprintf('  %-10s  %8.1f  %6.3f  %8.3f  %8.3f\n', components{ii}, M_gmol(ii), acentric(ii), 4.0, tau_sr(ii));
end

h_oil_base = [175, 204, 228, 327];
h_gas_base = [120, 100, 80, 60, 40, 20, 0];
measured_depths = [0; 175; 204; 228; 327];

results = struct();

%% ========================================================================
%  MODELS 1-4: Isothermal, Haase, Kempers, Firoozabadi (tau=4)
%  ========================================================================
model_ids   = [0, 1, 2, 3];
model_names = {'Isothermal', 'Haase', 'Kempers', 'Firoozabadi'};

for m = 1:4
    mid = model_ids(m);
    fprintf('\n=== Running %s ===\n', model_names{m});

    if mid == 3
        [GOC_h, ~, ~, GOC_T] = detect_sgoc(mid, comp_ref, press_ref, temp_ref, h_ref, dTdh, [0 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau);
    else
        [GOC_h, ~, ~, GOC_T] = detect_sgoc(mid, comp_ref, press_ref, temp_ref, h_ref, dTdh, [0 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    end

    h_oil = [GOC_h, GOC_h+2, h_oil_base];
    h_oil = unique(h_oil);
    n_oil = length(h_oil);
    oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1);
    oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil, n);

    for i = 1:n_oil
        try
            switch mid
                case 0, [ch,Ph,Th,~,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, 0, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 1, [ch,Ph,Th,~,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 2, [ch,Ph,Th,~,~] = main_kempers(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 3, [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
            end
            oil_P(i) = Ph; oil_comp(i,:) = ch';
            try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); oil_Pbub(i) = Pb; catch, end
            if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
                try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); oil_Pbub(i) = Pb; catch, end
            end
            if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
            try, oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, oil_rho(i) = NaN; end
        catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
    end

    [~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter);

    h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
    n_gas = length(h_gas);
    gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1);
    gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas, n);

    for i = 1:n_gas
        try
            switch mid
                case 0, [ch,Ph,Th,~,~] = main_hasse(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, 0, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 1, [ch,Ph,Th,~,~] = main_hasse(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 2, [ch,Ph,Th,~,~] = main_kempers(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
                case 3, [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
            end
            gas_P(i) = Ph; gas_comp(i,:) = ch';
            try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); gas_Pdew(i) = Pd; catch, end
            if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
                try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); gas_Pdew(i) = Pd; catch, end
            end
            if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
            try, gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, gas_rho(i) = NaN; end
        catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
    end

    all_h    = [flipud(h_gas'); h_oil'];
    all_comp = [flipud(gas_comp); oil_comp];

    results(m).name = model_names{m};
    results(m).GOC  = GOC_h;
    results(m).h    = all_h;
    results(m).P    = [flipud(gas_P); oil_P];
    results(m).Psat = [flipud(gas_Pdew); oil_Pbub];
    results(m).rho  = [flipud(gas_rho); oil_rho];
    results(m).comp = [flipud(gas_comp); oil_comp];
    results(m).C1   = all_comp(:, idx_C1) * 100;
    results(m).C7p  = sum(all_comp(:, idx_C7p), 2) * 100;
end

%% ========================================================================
%  MODEL 5: FIROOZABADI (SR param-free)
%  ========================================================================
fprintf('\n=== Running Firoozabadi (SR param-free) ===\n');

[GOC_h, ~, ~, GOC_T] = detect_sgoc(3, comp_ref, press_ref, temp_ref, h_ref, dTdh, [0 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau_sr(1));
fprintf('GOC: %.1f m\n', GOC_h);

h_oil = unique([GOC_h, GOC_h+2, h_oil_base]);
n_oil = length(h_oil);
oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1);
oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil, n);

for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_sr, vt_method, vt_params);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        try, oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, oil_rho(i) = NaN; end
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end

[~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter);

h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
n_gas = length(h_gas);
gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1);
gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas, n);

for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_sr, vt_method, vt_params);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        try, gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, gas_rho(i) = NaN; end
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end

all_h    = [flipud(h_gas'); h_oil'];
all_comp = [flipud(gas_comp); oil_comp];

results(5).name = 'Firoozabadi (SR)';
results(5).GOC  = GOC_h;
results(5).h    = all_h;
results(5).P    = [flipud(gas_P); oil_P];
results(5).Psat = [flipud(gas_Pdew); oil_Pbub];
results(5).rho  = [flipud(gas_rho); oil_rho];
results(5).comp = [flipud(gas_comp); oil_comp];
results(5).C1   = all_comp(:, idx_C1) * 100;
results(5).C7p  = sum(all_comp(:, idx_C7p), 2) * 100;

%% ========================================================================
%  EXPERIMENTAL DATA (Pedersen 2015, Table 1)
%  ========================================================================
exp_h    = [0; 175; 204; 228; 327];
exp_P    = [279; 284; 286; 287; 293];
exp_Psat = [270; 272; 267; 265; 242];
exp_C1   = [75.66; 50.04; 49.88; 48.89; 45.66];
exp_C7p  = [1.29+1.06+0.54+1.57; 3.50+3.75+2.28+15.88; 3.49+3.73+2.31+16.11; ...
            3.62+3.85+2.36+16.70; 4.01+4.28+2.58+17.66];

n_models = 5;
GOC_field = 140;

%% ========================================================================
%  PLOTS (2x2, clean)
%  ========================================================================
colors = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], ...
          [0.49 0.18 0.56], [0.00 0.75 0.75]};
lstyles = {'-', '--', ':', '-.', '-'};
lwidths = [1.5, 1.5, 1.5, 1.5, 2.0];

idx_plot = exp_h ~= h_ref;

fig = figure('Position', [50 50 1400 900], 'Color', 'w');

% ---- (a) C1 mol% ----
subplot(2,2,1); hold on;
for m = 1:n_models
    plot(results(m).C1, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C1(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse', 'Box', 'off'); grid on;
xlabel('C_1 (mol%)'); ylabel('Depth (m)');
text(0.03, 0.97, '(a)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

% ---- (b) C7+ mol% ----
subplot(2,2,2); hold on;
for m = 1:n_models
    plot(results(m).C7p, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C7p(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse', 'Box', 'off'); grid on;
xlabel('C_{7+} (mol%)'); ylabel('Depth (m)');
text(0.03, 0.97, '(b)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
legend([{results.name}, 'Experimental'], 'Location', 'best', 'FontSize', 6);

% ---- (c) Psat ----
subplot(2,2,3); hold on;
for m = 1:n_models
    plot(results(m).P/1e5, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m), 'HandleVisibility', 'off');
    v = ~isnan(results(m).Psat) & results(m).Psat > 0;
    plot(results(m).Psat(v)/1e5, results(m).h(v), 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_P(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(exp_Psat(idx_plot), exp_h(idx_plot), 'k^', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse', 'Box', 'off'); grid on;
xlabel('Pressure (bar)'); ylabel('Depth (m)');
text(0.03, 0.97, '(c)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
legend([{results.name}, 'P_{res}', 'P_{sat}'], 'Location', 'best', 'FontSize', 6);

% ---- (d) Density ----
subplot(2,2,4); hold on;
for m = 1:n_models
    v = ~isnan(results(m).rho) & results(m).rho > 0;
    plot(results(m).rho(v), results(m).h(v), 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse', 'Box', 'off'); grid on;
xlabel('Density (kg/m^3)'); ylabel('Depth (m)');
text(0.03, 0.97, '(d)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
legend({results.name}, 'Location', 'best', 'FontSize', 6);

sgtitle('Reservoir 1 (Pedersen 2015): 5-Model Comparison', 'FontSize', 13);

%% ========================================================================
%  GOC + RMSE TABLE
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  GOC COMPARISON (exp: ~140 m)\n');
fprintf('==========================================================\n');
for m = 1:n_models
    fprintf('  %-25s  GOC = %.1f m\n', results(m).name, results(m).GOC);
end

fprintf('\n==========================================================\n');
fprintf('  RMSE: OIL ZONE ONLY (excl. reference)\n');
fprintf('==========================================================\n');
fprintf('  %-25s  %8s  %8s  %9s\n', 'Model', 'RMSE(C1)', 'RMSE(C7+)', 'RMSE(all)');
fprintf('  %s\n', repmat('-', 1, 55));
is_gas = exp_h < GOC_field;
for m = 1:n_models
    err_C1 = []; err_C7p = []; err_all = [];
    for k = 1:length(exp_h)
        if is_gas(k), continue; end
        if exp_h(k) == h_ref, continue; end
        [~, idx] = min(abs(results(m).h - exp_h(k)));
        if abs(results(m).h(idx) - exp_h(k)) < 5
            err_C1  = [err_C1;  results(m).C1(idx) - exp_C1(k)];
            err_C7p = [err_C7p; results(m).C7p(idx) - exp_C7p(k)];
            err_all = [err_all; results(m).C1(idx) - exp_C1(k)];
            err_all = [err_all; results(m).C7p(idx) - exp_C7p(k)];
        end
    end
    if ~isempty(err_C1)
        fprintf('  %-25s  %7.2f %%  %7.2f %%  %9.2f %%\n', results(m).name, ...
            sqrt(mean(err_C1.^2)), sqrt(mean(err_C7p.^2)), sqrt(mean(err_all.^2)));
    end
end