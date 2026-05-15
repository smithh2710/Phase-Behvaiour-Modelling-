clear; clc; close all;

%% ========================================================================
%  PEDERSEN 2024 (SPE-223251) - HAASE REPLICATION + FIROOZABADI tau=4
%  SRK EOS + Peneloux volume shift, 14 components
%  Reference fluid at 2590 m
%  ========================================================================

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6', ...
              'C7-C10','C11-C20','C21-C30','C31-C80'};

n = 14;

comp_ref = [0.381; 0.654; 50.336; 7.838; 6.692; 1.019; 3.226; 1.148; 1.586; 1.941; ...
            11.911; 8.801; 2.839; 1.628] / 100;
comp_ref = comp_ref / sum(comp_ref);

M_gmol = [28.01; 44.01; 16.04; 30.07; 44.10; 58.12; 58.12; 72.15; 72.15; 86.18; ...
          110.24; 198.93; 341.69; 562.32];

Tc = [-146.95; 31.05; -82.55; 32.25; 96.65; 134.95; 152.05; 187.25; 196.45; 234.25; ...
      268.91; 388.66; 523.13; 719.41] + 273.15;

Pc = [33.94; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; 29.69; ...
      29.09; 19.59; 15.83; 14.66] * 1e5;

acentric = [0.0400; 0.2250; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; ...
            0.2510; 0.2960; 0.5135; 0.7671; 1.0736; 1.2903];

c_peneloux_cm3 = [0.920; 3.028; 0.630; 2.630; 5.060; 7.290; 7.860; 10.930; ...
                  12.180; 17.980; 9.707; 25.694; -1.051; -86.619];

R = 8.3144598;
Zc_est = 0.29056 - 0.08775 * acentric;
Vc = Zc_est .* R .* Tc ./ Pc;
Zc = Pc .* Vc ./ (R .* Tc);

BIP = zeros(n);
BIP(2,1) = -0.0315;
BIP(3,1) = 0.0278; BIP(4,1) = 0.0407; BIP(5,1) = 0.0763;
BIP(6,1) = 0.0944; BIP(7,1) = 0.0700; BIP(8,1) = 0.0867;
BIP(9,1) = 0.0878; BIP(10,1) = 0.0800;
for i = 11:14, BIP(i,1) = 0.0800; end
for i = 3:10,  BIP(i,2) = 0.1200; end
for i = 11:14, BIP(i,2) = 0.1000; end
BIP = BIP + BIP';

%% Cp coefficients (Poling for defined; Kesler-Lee for pseudo)
Cp_coeffs = [
    31.15,  -0.014,    2.68e-5,  -1.17e-8;   % N2
    19.79,   0.073,   -5.60e-5,   1.72e-8;   % CO2
    19.25,   0.052,    1.20e-5,  -1.13e-8;   % C1
     5.41,   0.178,   -6.94e-5,   8.71e-9;   % C2
    -4.22,   0.306,   -1.59e-4,   3.21e-8;   % C3
    -1.39,   0.385,   -1.85e-4,   2.90e-8;   % iC4
     9.49,   0.331,   -1.11e-4,  -2.82e-9;   % nC4
     9.52,   0.507,   -2.73e-4,   5.72e-8;   % iC5
    -3.63,   0.487,   -2.58e-4,   5.30e-8;   % nC5
    -4.41,   0.582,   -3.12e-4,   6.49e-8;   % C6
     4.50,   0.620,   -2.20e-4,   0.000;     % C7-C10
     7.50,   1.120,   -4.50e-4,   0.000;     % C11-C20
    12.00,   1.920,   -7.70e-4,   0.000;     % C21-C30
    20.00,   3.160,   -1.27e-3,   0.000];    % C31-C80

%% H_ig_ref: Pedersen 2024 Table B.1 TUNED values
%  Paper gives J/g (specific). Code needs J/mol (molar).
%  Conversion: H_ig_ref [J/mol] = H_ig [J/g] * M [g/mol]
H_ig_Jg = [0.0; 1231.4; 0.0; 304.3; 534.1; 620.3; 613.9; 461.5; 399.9; 393.3; ...
           147.1; 873.8; 883.2; 84.4];
H_ig_ref_tuned = H_ig_Jg .* M_gmol;

% Untuned: Pedersen 2006 correlation H_ig_ref/R = -1342 + 83.67*M
H_ig_ref_untuned = R * (-1342 + 83.67 * M_gmol);

%% Reservoir conditions
h_ref     = 2590;
press_ref = 284.5e5;
temp_ref  = 94.9 + 273.15;
dTdh      = 0.027;

tol = 1e-10;
maxiter = 500;
tau = 4.0;

%% VT setup: SRK + Peneloux (custom shifts)
vt_method = 12;
vt_params = struct();
vt_params.c_custom = c_peneloux_cm3;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

%% Component indices
idx_C1  = 3;
idx_C7p = 11:14;

%% ========================================================================
%  EXPERIMENTAL DATA: All 15 samples from Table A.1 (SPE-223251)
%  ========================================================================
exp_h = [2431; 2559; 2590; 2618; 2635; 2636; 2658; 2660; 2697; 2712; 2728; ...
         2746; 2758; 2758; 2769];

exp_C1 = [75.747; 75.228; 49.946; 49.510; 49.387; 49.884; 48.890; 47.936; 47.234; 47.213; 47.096; ...
          45.028; 45.323; 45.322; 44.882];

% C7+ = 100 - sum(N2 .. C6), computed from individual components in Table A.1
exp_light = [95.653; 95.463; 74.713; 74.276; 73.529; 74.360; 73.469; 72.636; 71.632; 72.069; 71.992; ...
             70.551; 70.933; 70.872; 70.348];
exp_C7p = 100 - exp_light;

exp_Pres = [278.3; 282.4; 284.2; 285.8; 286.6; 286.7; 288.1; 288.2; 290.4; 291.3; 292.4; ...
            291.6; 292.4; 292.4; 293.1];
exp_Psat = [278.3; 282.4; 284.1; 280.0; 268.5; 267.4; 265.0; 259.2; 254.2; 255.3; 253.8; ...
            241.2; 241.6; 240.6; 239.1];
exp_GOR  = [3875; 3061; 278.9; 276.4; NaN; 280.9; 264.5; NaN; 235.5; 229.3; 228.7; ...
            220.7; 227.8; 221.7; 213.4];

is_gas = abs(exp_Pres - exp_Psat) < 2;

% Individual C2-C6 from Table A.1 (mol%)
exp_C2 = [7.629; 8.001; 7.877; 7.828; 7.747; 7.775; 7.742; 7.568; 7.511; 7.588; 7.634; ...
          7.469; 7.541; 7.548; 7.475];
exp_C3 = [5.338; 5.547; 6.744; 6.746; 6.688; 6.745; 6.763; 6.799; 6.726; 6.824; 6.871; ...
          7.021; 7.073; 7.059; 7.048];
exp_iC4 = [0.735; 0.737; 1.031; 1.031; 1.008; 1.035; 1.044; 1.050; 1.038; 1.063; 1.056; ...
           1.114; 1.119; 1.111; 1.117];
exp_nC4 = [2.023; 2.092; 3.265; 3.253; 3.181; 3.206; 3.241; 3.317; 3.276; 3.352; 3.332; ...
           3.475; 3.492; 3.453; 3.480];
exp_iC5 = [0.669; 0.647; 1.190; 1.182; 1.082; 1.162; 1.181; 1.158; 1.142; 1.169; 1.160; ...
           1.306; 1.300; 1.275; 1.298];
exp_nC5 = [0.837; 0.826; 1.644; 1.627; 1.483; 1.551; 1.581; 1.579; 1.557; 1.589; 1.577; ...
           1.743; 1.739; 1.696; 1.733];
exp_C6 = [0.911; 0.850; 2.015; 2.002; 1.882; 1.878; 1.934; 2.019; 1.996; 2.008; 1.997; ...
          2.187; 2.158; 2.212; 2.158];
exp_C2C6 = [exp_C2, exp_C3, exp_iC4, exp_nC4, exp_iC5, exp_nC5, exp_C6];
idx_C2C6 = [4, 5, 6, 7, 8, 9, 10];

% Split C7+ from Table A.1: C7-C10
exp_C7_C10 = [1.363+1.087+0.546+0.369; 1.374+1.176+0.557+0.378; ...
              4.006+3.969+2.303+1.822; 3.994+3.996+2.335+1.854; ...
              3.958+4.089+2.464+1.664; 3.719+3.840+2.358+1.854; ...
              3.851+3.970+2.408+1.940; 4.123+4.237+2.571+1.990; ...
              4.183+4.324+2.612+1.794; 4.096+4.207+2.605+1.996; ...
              4.069+4.215+2.593+2.016; 4.281+4.440+2.647+2.206; ...
              4.223+4.373+2.616+2.158; 4.190+4.260+2.589+2.184; ...
              4.229+4.387+2.630+2.197];

% C11-C19
exp_C11_C19 = [0.234+0.173+0.140+0.106+0.089+0.060+0.046+0.039+0.028; ...
               0.242+0.175+0.150+0.112+0.096+0.067+0.053+0.046+0.033; ...
               1.346+1.150+1.134+0.976+0.979+0.782+0.673+0.644+0.535; ...
               1.376+1.178+1.164+0.998+1.010+0.804+0.696+0.662+0.550; ...
               1.491+1.336+1.196+1.072+0.960+0.860+0.770+0.690+0.618; ...
               1.427+1.231+1.163+1.021+1.046+0.810+0.720+0.650+0.561; ...
               1.480+1.275+1.203+1.054+1.086+0.839+0.746+0.673+0.583; ...
               1.549+1.336+1.263+1.121+1.073+0.873+0.733+0.731+0.629; ...
               1.607+1.440+1.291+1.156+1.036+0.928+0.832+0.745+0.668; ...
               1.572+1.352+1.301+1.155+1.097+0.900+0.757+0.753+0.653; ...
               1.601+1.355+1.300+1.153+1.104+0.913+0.763+0.763+0.668; ...
               1.626+1.407+1.283+1.190+1.202+0.934+0.822+0.789+0.664; ...
               1.619+1.373+1.275+1.168+1.191+0.932+0.836+0.781+0.638; ...
               1.616+1.403+1.282+1.190+1.210+0.945+0.825+0.789+0.667; ...
               1.629+1.417+1.296+1.206+1.262+0.952+0.837+0.814+0.678];

exp_C20p = [0.064; 0.078; 4.968; 5.107; 5.303; 5.240; 5.423; 5.135; 5.751; 5.487; 5.495; ...
            5.958; 5.884; 5.978; 6.118];

exp_C11p = exp_C11_C19 + exp_C20p;

%% ========================================================================
%  RUN 5 MODELS
%  ========================================================================

Tc_mix = sum(comp_ref .* Tc);
Tr_mix = temp_ref / Tc_mix;

% Model 4: SR parameter-free: tau = Tr^(1/4) + ln(M)*omega^2 - omega
tau_sr = Tr_mix^0.25 + log(M_gmol) .* acentric.^2 - acentric;
tau_sr = max(min(tau_sr, 10.0), 0.9);

% Model 5: Regressed 5-param: tau = a*Tr + b*ln(M) + c*omega + d*ln(M)*omega + e*Tr*omega
p_reg = [-2.5542, -0.2679, 1.1697, 1.8033, 0.4798];
tau_reg = p_reg(1)*Tr_mix + p_reg(2)*log(M_gmol) + p_reg(3)*acentric + ...
          p_reg(4)*log(M_gmol).*acentric + p_reg(5)*Tr_mix*acentric;
tau_reg = max(min(tau_reg, 10.0), 0.9);

fprintf('=== Tau correlation summary ===\n');
fprintf('  Tr_mix = %.4f,  Tr^(1/4) = %.4f\n', Tr_mix, Tr_mix^0.25);
fprintf('  %-10s  %8s  %6s  %8s  %8s  %8s\n', 'Component', 'M', 'omega', 'tau=4', 'tau_SR', 'tau_reg');
fprintf('  %s\n', repmat('-', 1, 58));
for i = 1:n
    fprintf('  %-10s  %8.1f  %6.3f  %8.3f  %8.3f  %8.3f\n', components{i}, M_gmol(i), acentric(i), 4.0, tau_sr(i), tau_reg(i));
end

n_models = 5;
model_names = {'Haase (Ped. 2006 H_{ig})', 'Haase (tuned H_{ig})', ...
               'Firoozabadi (\tau=4)', 'Firoozabadi (SR param-free)', ...
               'Firoozabadi (regressed \tau)'};

% Per-model config
model_cfg = struct();
model_cfg(1).type = 1;  model_cfg(1).dT = dTdh;  model_cfg(1).H = H_ig_ref_untuned; model_cfg(1).tau = [];
model_cfg(2).type = 1;  model_cfg(2).dT = dTdh;  model_cfg(2).H = H_ig_ref_tuned;   model_cfg(2).tau = [];
model_cfg(3).type = 3;  model_cfg(3).dT = dTdh;  model_cfg(3).H = H_ig_ref_tuned;   model_cfg(3).tau = tau;
model_cfg(4).type = 3;  model_cfg(4).dT = dTdh;  model_cfg(4).H = H_ig_ref_tuned;   model_cfg(4).tau = tau_sr;
model_cfg(5).type = 3;  model_cfg(5).dT = dTdh;  model_cfg(5).H = H_ig_ref_tuned;   model_cfg(5).tau = tau_reg;

h_oil_base = unique([2590, 2600, 2618, 2635, 2636, 2658, 2660, 2697, 2712, 2728, 2746, 2758, 2769]);
h_gas_base = [2559, 2540, 2500, 2460, 2431];

results = struct();

for m = 1:n_models
    mid  = model_cfg(m).type;
    dT_m = model_cfg(m).dT;
    H_m  = model_cfg(m).H;
    tau_m = model_cfg(m).tau;
    fprintf('\n=== %s (type=%d) ===\n', model_names{m}, mid);

    if mid == 3
        [GOC_h, ~, ~, GOC_T] = detect_sgoc(mid, comp_ref, press_ref, temp_ref, h_ref, dTdh, ...
            [2500 2650], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, vt_method, vt_params, tau_m);
    else
        [GOC_h, ~, ~, GOC_T] = detect_sgoc(mid, comp_ref, press_ref, temp_ref, h_ref, dT_m, ...
            [2500 2650], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, vt_method, vt_params);
    end
    fprintf('  GOC = %.2f m\n', GOC_h);

    h_oil = [GOC_h, GOC_h+2, h_oil_base];
    n_oil = length(h_oil);
    oil_P = zeros(n_oil,1); oil_Pbub = zeros(n_oil,1);
    oil_rho = zeros(n_oil,1); oil_comp = zeros(n_oil, n);
    oil_GOR = NaN(n_oil,1);

    for i = 1:n_oil
        if mid == 3
            [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, tau_m, vt_method, vt_params);
        else
            [ch,Ph,Th,~,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dT_m, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, vt_method, vt_params);
        end
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try
            [Pb_ext, ~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
            oil_Pbub(i) = Pb_ext;
        catch, oil_Pbub(i) = NaN;
        end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try
                [Pb_ext, ~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
                oil_Pbub(i) = Pb_ext;
            catch, oil_Pbub(i) = NaN;
            end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph)
            oil_Pbub(i) = NaN;
        end
        try
            oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
        catch, oil_rho(i) = NaN;
        end
        try
            oil_GOR(i) = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, ...
                'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'SRK');
        catch, oil_GOR(i) = NaN;
        end
    end

    [~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, ...
        Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

    h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
    n_gas = length(h_gas);
    gas_P = zeros(n_gas,1); gas_Pdew = zeros(n_gas,1);
    gas_rho = zeros(n_gas,1); gas_comp = zeros(n_gas, n);
    gas_GOR = NaN(n_gas,1);

    for i = 1:n_gas
        if mid == 3
            [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, tau_m, vt_method, vt_params);
        else
            [ch,Ph,Th,~,~] = main_hasse(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dT_m, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_m, vt_method, vt_params);
        end
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try
            [Pd_ext, ~] = pressdew_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
            gas_Pdew(i) = Pd_ext;
        catch, gas_Pdew(i) = NaN;
        end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try
                [Pd_ext, ~] = pressdew_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
                gas_Pdew(i) = Pd_ext;
            catch, gas_Pdew(i) = NaN;
            end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph)
            gas_Pdew(i) = NaN;
        end
        try
            gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
        catch, gas_rho(i) = NaN;
        end
        try
            gas_GOR(i) = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, ...
                'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'SRK');
        catch, gas_GOR(i) = NaN;
        end
    end

    all_h    = [flipud(h_gas'); h_oil'];
    all_comp = [flipud(gas_comp); oil_comp];
    all_P    = [flipud(gas_P); oil_P];
    all_Psat = [flipud(gas_Pdew); oil_Pbub];
    all_rho  = [flipud(gas_rho); oil_rho];
    all_GOR  = [flipud(gas_GOR); oil_GOR];

    results(m).name = model_names{m};
    results(m).GOC  = GOC_h;
    results(m).h    = all_h;
    results(m).P    = all_P;
    results(m).Psat = all_Psat;
    results(m).rho  = all_rho;
    results(m).C1   = all_comp(:, idx_C1) * 100;
    results(m).C7p  = sum(all_comp(:, idx_C7p), 2) * 100;
    results(m).comp = all_comp;
    results(m).GOR  = all_GOR;
end

%% ========================================================================
%  PLOTS (matching Pedersen 2024 Figure 3 format)
%  ========================================================================

colors = {[0.47 0.67 0.19], ...   % green   - Haase untuned
          [0.85 0.33 0.10], ...   % red     - Haase tuned
          [0.49 0.18 0.56], ...   % purple  - Firoozabadi tau=4
          [0.93 0.69 0.13], ...   % gold    - SR param-free
          [0.00 0.75 0.75]};      % cyan    - Regressed

lstyles = {':', '--', '-.', '-', '-'};
lwidths = [1.5, 1.5, 1.5, 1.8, 2.0];
figure('Position', [50 50 1400 900]);

% Exclude reference depth (2590 m) from experimental plots
idx_plot = exp_h ~= 2590;

% (a) C1 mol%
subplot(2,2,1);
hold on;
for m = 1:n_models
    plot(results(m).C1, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C1(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('C_1 (mol%)'); ylabel('Depth (m)');
legend([model_names, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(a) C_1 mol% vs. Depth');

% (b) C7+ mol%
subplot(2,2,2);
hold on;
for m = 1:n_models
    plot(results(m).C7p, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C7p(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('C_{7+} (mol%)'); ylabel('Depth (m)');
legend([model_names, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(b) C_{7+} mol% vs. Depth');

% (c) Pressure + Saturation Pressure
subplot(2,2,3);
hold on;
for m = 1:n_models
    plot(results(m).Psat/1e5, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(results(1).P/1e5, results(1).h, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(exp_Pres(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
plot(exp_Psat(idx_plot), exp_h(idx_plot), 'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('Pressure (bar)'); ylabel('Depth (m)');
legend([model_names, 'P_{res}', 'P_{sat}'], 'Location', 'best', 'FontSize', 6);
title('(c) P^{res} and P^{sat} vs. Depth');

% (d) GOR
subplot(2,2,4);
hold on;
for m = 1:n_models
    valid = ~isnan(results(m).GOR) & results(m).GOR > 0;
    if any(valid)
        plot(results(m).GOR(valid), results(m).h(valid), 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
    end
end
valid_gor = ~isnan(exp_GOR) & idx_plot;
plot(exp_GOR(valid_gor), exp_h(valid_gor), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('GOR (Sm^3/Sm^3)'); ylabel('Depth (m)');
legend([model_names, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(d) GOR vs. Depth');

sgtitle('Pedersen 2024 (SPE-223251): 5-Model Comparison', 'FontSize', 13);

%% ========================================================================
%  QUANTITATIVE COMPARISON
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  GOC COMPARISON\n');
fprintf('==========================================================\n');
fprintf('  Experimental GOC range:       2559-2590 m\n');
fprintf('  Pedersen 2024 simulated GOC:  2571 m\n');
for m = 1:n_models
    fprintf('  %-30s  GOC = %.1f m\n', results(m).name, results(m).GOC);
end

fprintf('\n==========================================================\n');
fprintf('  RMSE: OIL ZONE ONLY (%d samples)\n', sum(~is_gas));
fprintf('==========================================================\n');
fprintf('  %-30s  %8s  %8s  %9s  %8s  %9s  %9s\n', 'Model', 'RMSE(C1)', 'RMSE(C7+)', 'RMSE(C7-10)', 'RMSE(C11+)', 'RMSE(C2-C6)', 'RMSE(all)');
fprintf('  %s\n', repmat('-', 1, 105));
for m = 1:n_models
    err_C1 = []; err_C7p = []; err_C710 = []; err_C11p = []; err_C2C6_all = []; err_all = [];
    for k = 1:length(exp_h)
        if is_gas(k), continue; end
        [~, idx] = min(abs(results(m).h - exp_h(k)));
        if abs(results(m).h(idx) - exp_h(k)) < 15
            err_C1   = [err_C1;   results(m).C1(idx) - exp_C1(k)];
            err_C7p  = [err_C7p;  results(m).C7p(idx) - exp_C7p(k)];
            err_C710 = [err_C710; results(m).comp(idx,11)*100 - exp_C7_C10(k)];
            err_C11p = [err_C11p; sum(results(m).comp(idx,12:14))*100 - exp_C11p(k)];
            for j = 1:length(idx_C2C6)
                err_C2C6_all = [err_C2C6_all; results(m).comp(idx,idx_C2C6(j))*100 - exp_C2C6(k,j)];
            end
            % Overall: all individual component errors at this depth
            err_all = [err_all; results(m).C1(idx) - exp_C1(k)];
            err_all = [err_all; results(m).comp(idx,11)*100 - exp_C7_C10(k)];
            err_all = [err_all; sum(results(m).comp(idx,12:14))*100 - exp_C11p(k)];
            for j = 1:length(idx_C2C6)
                err_all = [err_all; results(m).comp(idx,idx_C2C6(j))*100 - exp_C2C6(k,j)];
            end
        end
    end
    if ~isempty(err_C1)
        fprintf('  %-30s  %7.2f %%  %7.2f %%  %9.2f %%  %8.2f %%  %9.2f %%  %9.2f %%\n', results(m).name, ...
            sqrt(mean(err_C1.^2)), sqrt(mean(err_C7p.^2)), sqrt(mean(err_C710.^2)), ...
            sqrt(mean(err_C11p.^2)), sqrt(mean(err_C2C6_all.^2)), sqrt(mean(err_all.^2)));
    end
end

fprintf('\n==========================================================\n');
fprintf('  SAMPLE-BY-SAMPLE COMPARISON (oil zone)\n');
fprintf('==========================================================\n');
fprintf('  %-6s | %8s | ', 'Depth', 'C1_exp');
for m = 1:n_models, fprintf('%8s | ', results(m).name(1:min(8,end))); end
fprintf('\n  %s\n', repmat('-', 1, 55));
for k = 1:length(exp_h)
    if is_gas(k), continue; end
    fprintf('  %-6.0f | %8.2f | ', exp_h(k), exp_C1(k));
    for m = 1:n_models
        [~, idx] = min(abs(results(m).h - exp_h(k)));
        if abs(results(m).h(idx) - exp_h(k)) < 15
            fprintf('%8.2f | ', results(m).C1(idx));
        else
            fprintf('%8s | ', 'N/A');
        end
    end
    fprintf('\n');
end