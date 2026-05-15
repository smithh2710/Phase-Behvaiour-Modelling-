clear; clc; close all;

%% ========================================================================
%  COMMON DATA (Reservoir 2 - Pedersen & Hjermstad 2006)
%  ========================================================================
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
press_ref = 380.20e5;
temp_ref = 138.70 + 273.15;
dTdh_field = 0.0260;
R = 8.3144598;

Zc = 0.29056 - 0.08775 .* acentric;
for i = 2:length(Zc)
    if Zc(i) >= Zc(i-1)
        Zc(i) = Zc(i-1) - 0.008;
    end
end
Vc = Zc .* R .* Tc ./ Pc;

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

tol = 1e-10;
maxiter = 500;
measured_depths = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];

comp_gas_GOC_hasse = [0.566; 2.279; 66.572; 8.069; 5.434; 0.940; 2.248; 0.766; 1.163; 1.226; ...
                      1.809; 1.849; 1.192; 2.518; 1.206; 0.792; 0.512; 0.402; 0.247; 0.135; 0.054; 0.021];
comp_gas_GOC_hasse = comp_gas_GOC_hasse ./ sum(comp_gas_GOC_hasse);

results = struct();
idx_C1  = 3;
idx_C7p = 11:22;

% SR param-free tau: tau = Tr^(1/4) + ln(M)*omega^2 - omega
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

%% ========================================================================
%  MODEL 1: ISOTHERMAL (main_hasse with dTdh = 0)
%  ========================================================================
fprintf('\n=== Running Isothermal ===\n');
vm = 7; vp.Vc = Vc; vp.components = components;

[GOC_h, ~, ~, ~] = detect_sgoc(1, comp_ref, press_ref, temp_ref, h_ref, 0, ...
    [3580 3682], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
fprintf('GOC: %.1f m\n', GOC_h);

oil_m = measured_depths(measured_depths > GOC_h + 2);
h_oil = sort(unique([GOC_h; GOC_h+10; oil_m]));
gas_m = measured_depths(measured_depths < GOC_h - 2);
h_gas = sort(unique([GOC_h; GOC_h-5; GOC_h-10; gas_m; 3600]), 'descend');
n_oil = length(h_oil); n_gas = length(h_gas);

oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1); oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, 0, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end
oil_Pbub(1) = oil_P(1);
GOC_T = temp_ref;
[~, cg_ref] = pressbub_multicomp_ss(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1); gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_hasse(h_gas(i), GOC_h, cg_ref, oil_P(1), GOC_T, 0, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end
gas_Pdew(1) = oil_P(1); gas_P(1) = oil_P(1);

results(1).name = 'Isothermal'; results(1).GOC = GOC_h;
results(1).h = [flipud(h_gas); h_oil]; results(1).P = [flipud(gas_P); oil_P];
results(1).Psat = [flipud(gas_Pdew); oil_Pbub]; results(1).rho = [flipud(gas_rho); oil_rho];
ac = [flipud(gas_comp); oil_comp];
results(1).comp = ac;
results(1).C1 = ac(:,idx_C1)*100; results(1).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 2: HAASE (hardcoded gas GOC comp)
%  ========================================================================
fprintf('\n=== Running Haase ===\n');
vm = 7; vp = struct('c_custom', zeros(n,1), 'Vc', Vc, 'components', {components}, 'Zc', Zc);

[GOC_h, ~, ~, ~] = detect_sgoc(1, comp_ref, press_ref, temp_ref, h_ref, dTdh_field, ...
    [3600 3680], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
fprintf('GOC: %.1f m\n', GOC_h);
GOC_T = temp_ref + dTdh_field * (GOC_h - h_ref);

h_oil = [GOC_h; 3655; 3661.6; 3670; 3676.0; 3682.8];
h_gas = [GOC_h; 3644.3; 3642; 3638.2; 3630; 3615; 3600];
n_oil = length(h_oil); n_gas = length(h_gas);

oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1); oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end
P_GOC = oil_P(1);
oil_Pbub(1) = oil_P(1);

gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1); gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_hasse(h_gas(i), GOC_h, comp_gas_GOC_hasse, P_GOC, GOC_T, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end
gas_Pdew(1) = oil_P(1); gas_P(1) = oil_P(1);

results(2).name = 'Haase'; results(2).GOC = GOC_h;
results(2).h = [flipud(h_gas); h_oil]; results(2).P = [flipud(gas_P); oil_P];
results(2).Psat = [flipud(gas_Pdew); oil_Pbub]; results(2).rho = [flipud(gas_rho); oil_rho];
ac = [flipud(gas_comp); oil_comp];
results(2).comp = ac;
results(2).C1 = ac(:,idx_C1)*100; results(2).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 3: KEMPERS (flash at GOC, dynamic depths)
%  ========================================================================
fprintf('\n=== Running Kempers ===\n');
vm = 7; vp = struct('c_custom', zeros(n,1), 'Vc', Vc, 'components', {components}, 'Zc', Zc);

[GOC_h, ~, ~, ~] = detect_sgoc(2, comp_ref, press_ref, temp_ref, h_ref, dTdh_field, ...
    [3580 3682], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
fprintf('GOC: %.1f m\n', GOC_h);
GOC_T = temp_ref + dTdh_field * (GOC_h - h_ref);

oil_m = measured_depths(measured_depths > GOC_h + 2);
h_oil = sort(unique([GOC_h; GOC_h+5; GOC_h+10; oil_m]));
gas_m = measured_depths(measured_depths < GOC_h - 2);
h_gas = sort(unique([GOC_h; GOC_h-5; GOC_h-10; gas_m; 3600]), 'descend');
n_oil = length(h_oil); n_gas = length(h_gas);

oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1); oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_kempers(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end
oil_Pbub(1) = oil_P(1);
[~, cg_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1); gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_kempers(h_gas(i), GOC_h, cg_ref, oil_P(1), GOC_T, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end
gas_Pdew(1) = oil_P(1); gas_P(1) = oil_P(1);

results(3).name = 'Kempers'; results(3).GOC = GOC_h;
results(3).h = [flipud(h_gas); h_oil]; results(3).P = [flipud(gas_P); oil_P];
results(3).Psat = [flipud(gas_Pdew); oil_Pbub]; results(3).rho = [flipud(gas_rho); oil_rho];
ac = [flipud(gas_comp); oil_comp];
results(3).comp = ac;
results(3).C1 = ac(:,idx_C1)*100; results(3).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 4: FIROOZABADI (flash at GOC, dynamic depths, tau=4)
%  ========================================================================
fprintf('\n=== Running Firoozabadi ===\n');
tau = 4;
vm = 7; vp.Vc = Vc; vp.components = components;

[GOC_h, ~, ~, ~] = detect_sgoc(3, comp_ref, press_ref, temp_ref, h_ref, dTdh_field, ...
    [3580 3682], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);
fprintf('GOC: %.1f m\n', GOC_h);
GOC_T = temp_ref + dTdh_field * (GOC_h - h_ref);

oil_m = measured_depths(measured_depths > GOC_h + 2);
h_oil = sort(unique([GOC_h; GOC_h+10; GOC_h+15; GOC_h+20; oil_m]));
gas_m = measured_depths(measured_depths < GOC_h - 2);
h_gas = sort(unique([GOC_h; GOC_h-5; GOC_h-10; GOC_h-15; GOC_h-20; gas_m; 3600]), 'descend');
n_oil = length(h_oil); n_gas = length(h_gas);

oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1); oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vm, vp);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end
oil_Pbub(1) = oil_P(1);
[~, cg_ref] = pressbub_multicomp_ss(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1); gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, cg_ref, oil_P(1), GOC_T, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vm, vp);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end
gas_Pdew(1) = oil_P(1); gas_P(1) = oil_P(1);

results(4).name = 'Firoozabadi'; results(4).GOC = GOC_h;
results(4).h = [flipud(h_gas); h_oil]; results(4).P = [flipud(gas_P); oil_P];
results(4).Psat = [flipud(gas_Pdew); oil_Pbub]; results(4).rho = [flipud(gas_rho); oil_rho];
ac = [flipud(gas_comp); oil_comp];
results(4).comp = ac;
results(4).C1 = ac(:,idx_C1)*100; results(4).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 5: FIROOZABADI (SR param-free: Tr^(1/4) + ln(M)*omega^2 - omega)
%  ========================================================================
fprintf('\n=== Running Firoozabadi (SR param-free) ===\n');
vm = 7; vp.Vc = Vc; vp.components = components;

[GOC_h, ~, ~, ~] = detect_sgoc(3, comp_ref, press_ref, temp_ref, h_ref, dTdh_field, ...
    [3580 3682], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp, tau_sr(1));
fprintf('GOC: %.1f m\n', GOC_h);
GOC_T = temp_ref + dTdh_field * (GOC_h - h_ref);

oil_m = measured_depths(measured_depths > GOC_h + 2);
h_oil = sort(unique([GOC_h; GOC_h+10; GOC_h+15; GOC_h+20; oil_m]));
gas_m = measured_depths(measured_depths < GOC_h - 2);
h_gas = sort(unique([GOC_h; GOC_h-5; GOC_h-10; GOC_h-15; GOC_h-20; gas_m; 3600]), 'descend');
n_oil = length(h_oil); n_gas = length(h_gas);

oil_P = NaN(n_oil,1); oil_Pbub = NaN(n_oil,1); oil_rho = NaN(n_oil,1); oil_comp = zeros(n_oil,n);
for i = 1:n_oil
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_sr, vm, vp);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try, [Pb,~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        if isnan(oil_Pbub(i)) || oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph
            try, [Pb,~] = pressbub_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); oil_Pbub(i) = Pb; catch, end
        end
        if ~isnan(oil_Pbub(i)) && (oil_Pbub(i) <= 0 || oil_Pbub(i) > Ph), oil_Pbub(i) = NaN; end
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Oil error at h=%.1f: %s\n', h_oil(i), ME.message); end
end
oil_Pbub(1) = oil_P(1);
[~, cg_ref] = pressbub_multicomp_ss(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

gas_P = NaN(n_gas,1); gas_Pdew = NaN(n_gas,1); gas_rho = NaN(n_gas,1); gas_comp = zeros(n_gas,n);
for i = 1:n_gas
    try
        [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, cg_ref, oil_P(1), GOC_T, dTdh_field, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_sr, vm, vp);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try, [Pd,~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        if isnan(gas_Pdew(i)) || gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph
            try, [Pd,~] = pressdew_multicomp_newton(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components); gas_Pdew(i) = Pd; catch, end
        end
        if ~isnan(gas_Pdew(i)) && (gas_Pdew(i) <= 0 || gas_Pdew(i) > Ph), gas_Pdew(i) = NaN; end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vm, vp);
    catch ME, fprintf('  Gas error at h=%.1f: %s\n', h_gas(i), ME.message); end
end
gas_Pdew(1) = oil_P(1); gas_P(1) = oil_P(1);

results(5).name = 'Firoozabadi (SR)'; results(5).GOC = GOC_h;
results(5).h = [flipud(h_gas); h_oil]; results(5).P = [flipud(gas_P); oil_P];
results(5).Psat = [flipud(gas_Pdew); oil_Pbub]; results(5).rho = [flipud(gas_rho); oil_rho];
ac = [flipud(gas_comp); oil_comp];
results(5).comp = ac;
results(5).C1 = ac(:,idx_C1)*100; results(5).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  EXPERIMENTAL DATA
%  ========================================================================
exp_h    = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];
exp_P    = [377.8; 377.9; 378.2; 378.8; 379.6; 380.2];
exp_Psat = [375.5; 372.8; 364.2; 364.5; 360.1; 353.8];
exp_rho  = [367; 376; 574; 581; 595; 601];
exp_C1   = [68.861; 68.546; 56.142; 55.261; 54.253; 53.871];
exp_C7p  = [1.482+1.595+1.031+5.231; 1.519+1.610+1.048+5.080; ...
            2.474+2.583+1.695+13.78; 2.520+2.667+1.779+14.491; ...
            2.579+2.777+1.869+15.282; 2.630+2.823+1.897+15.783];

%% ========================================================================
%  PLOTS (2x2 layout)
%  ========================================================================
n_models = 5;
colors = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], ...
          [0.49 0.18 0.56], [0.00 0.75 0.75]};
lstyles = {'-', '--', ':', '-.', '-'};
lwidths = [1.5, 1.5, 1.5, 1.5, 2.0];
idx_plot = exp_h ~= h_ref;
GOC_field = 3647;

figure('Position', [50 50 1400 900], 'Color', 'w');

subplot(2,2,1); hold on;
for m = 1:n_models
    plot(results(m).C1, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C1(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse'); grid on; ylim([3600 3700]);
xlabel('C_1 (mol%)'); ylabel('Depth (m)');
legend([{results.name}, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(a) C_1 mol% vs. Depth');

subplot(2,2,2); hold on;
for m = 1:n_models
    plot(results(m).C7p, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C7p(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse'); grid on; ylim([3600 3700]);
xlabel('C_{7+} (mol%)'); ylabel('Depth (m)');
legend([{results.name}, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(b) C_{7+} mol% vs. Depth');

subplot(2,2,3); hold on;
for m = 1:n_models
    plot(results(m).Psat/1e5, results(m).h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(results(1).P/1e5, results(1).h, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(exp_P(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
plot(exp_Psat(idx_plot), exp_h(idx_plot), 'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse'); grid on; ylim([3600 3700]);
xlabel('Pressure (bar)'); ylabel('Depth (m)');
legend([{results.name}, 'P_{res}', 'P_{sat}'], 'Location', 'best', 'FontSize', 6);
title('(c) P^{res} and P^{sat} vs. Depth');

subplot(2,2,4); hold on;
for m = 1:n_models
    v = ~isnan(results(m).rho) & results(m).rho > 0;
    plot(results(m).rho(v), results(m).h(v), 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_rho(idx_plot), exp_h(idx_plot), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
yline(GOC_field, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(gca, 'YDir', 'reverse'); grid on; ylim([3600 3700]);
xlabel('Density (kg/m^3)'); ylabel('Depth (m)');
legend([{results.name}, 'Experimental'], 'Location', 'best', 'FontSize', 6);
title('(d) Density vs. Depth');

sgtitle('Pedersen 2006 (SPE-101275): 5-Model Comparison', 'FontSize', 13);

%% ========================================================================
%  RMSE TABLE
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  GOC COMPARISON (exp: ~3647 m)\n');
fprintf('==========================================================\n');
for m = 1:n_models
    fprintf('  %-30s  GOC = %.1f m\n', results(m).name, results(m).GOC);
end

fprintf('\n==========================================================\n');
fprintf('  RMSE: OIL ZONE ONLY (excl. reference)\n');
fprintf('==========================================================\n');

% Full experimental compositions from SPE-101275 Table 1 [mol%]
% Rows: depths [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8]
% Cols: N2, CO2, C1, C2, C3, iC4, nC4, iC5, nC5, C6, C7, C8, C9, C10+
exp_comp = [
    0.431, 2.752, 68.861, 8.427, 5.198, 0.847, 1.885, 0.587, 0.752, 0.921, 1.482, 1.595, 1.031, 5.231;
    0.295, 2.834, 68.546, 8.341, 5.212, 0.892, 2.100, 0.675, 0.866, 0.981, 1.519, 1.610, 1.048, 5.080;
    0.358, 2.332, 56.142, 8.094, 5.535, 1.001, 2.439, 0.879, 1.184, 1.504, 2.474, 2.583, 1.695, 13.78;
    0.331, 2.455, 55.261, 8.025, 5.481, 0.995, 2.433, 0.877, 1.182, 1.504, 2.520, 2.667, 1.779, 14.491;
    0.337, 2.363, 54.253, 7.961, 5.494, 1.000, 2.454, 0.889, 1.202, 1.539, 2.579, 2.777, 1.869, 15.282;
    0.395, 2.060, 53.871, 7.589, 5.575, 1.009, 2.514, 0.900, 1.396, 1.557, 2.630, 2.823, 1.897, 15.783];

% Model component indices mapping to exp_comp columns:
%  exp col 1 (N2)   -> model comp 1
%  exp col 2 (CO2)  -> model comp 2
%  exp col 3 (C1)   -> model comp 3
%  exp col 4-10 (C2-C6) -> model comp 4-10
%  exp col 11 (C7)  -> model comp 11
%  exp col 12 (C8)  -> model comp 12
%  exp col 13 (C9)  -> model comp 13
%  exp col 14 (C10+)-> sum(model comp 14:22)

is_gas = abs(exp_P - exp_Psat) < 3;

fprintf('  %-25s  %7s  %7s  %7s  %9s  %8s  %8s\n', ...
    'Model', 'RMSE(N2)', 'RMSE(CO2)', 'RMSE(C1)', 'RMSE(C2-C6)', 'RMSE(C7+)', 'RMSE(all)');
fprintf('  %s\n', repmat('-', 1, 85));

for m = 1:n_models
    err_N2 = []; err_CO2 = []; err_C1 = []; err_C2C6 = []; err_C7p = []; err_all = [];
    for k = 1:length(exp_h)
        if is_gas(k), continue; end
        if exp_h(k) == h_ref, continue; end
        [~, idx] = min(abs(results(m).h - exp_h(k)));
        if abs(results(m).h(idx) - exp_h(k)) > 5, continue; end

        mc = results(m).comp(idx,:) * 100;

        % N2
        e = mc(1) - exp_comp(k,1);
        err_N2 = [err_N2; e]; err_all = [err_all; e];

        % CO2
        e = mc(2) - exp_comp(k,2);
        err_CO2 = [err_CO2; e]; err_all = [err_all; e];

        % C1
        e = mc(3) - exp_comp(k,3);
        err_C1 = [err_C1; e]; err_all = [err_all; e];

        % C2-C6 (individual components 4-10 vs exp cols 4-10)
        for j = 4:10
            e = mc(j) - exp_comp(k,j);
            err_C2C6 = [err_C2C6; e]; err_all = [err_all; e];
        end

        % C7+ (model sum 11:22 vs exp C7+C8+C9+C10+)
        mod_C7p = sum(mc(11:22));
        exp_C7p_k = sum(exp_comp(k,11:14));
        e = mod_C7p - exp_C7p_k;
        err_C7p = [err_C7p; e]; err_all = [err_all; e];
    end
    if ~isempty(err_C1)
        fprintf('  %-25s  %7.2f%%  %8.2f%%  %7.2f%%  %10.2f%%  %8.2f%%  %8.2f%%\n', results(m).name, ...
            sqrt(mean(err_N2.^2)), sqrt(mean(err_CO2.^2)), sqrt(mean(err_C1.^2)), ...
            sqrt(mean(err_C2C6.^2)), sqrt(mean(err_C7p.^2)), sqrt(mean(err_all.^2)));
    end
end