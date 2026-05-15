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

%% ========================================================================
%  MODEL 1: ISOTHERMAL
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
results(1).C1 = ac(:,idx_C1)*100; results(1).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 2: HAASE
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
results(2).C1 = ac(:,idx_C1)*100; results(2).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 3: KEMPERS
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
results(3).C1 = ac(:,idx_C1)*100; results(3).C7p = sum(ac(:,idx_C7p),2)*100;

%% ========================================================================
%  MODEL 4: FIROOZABADI
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
results(4).C1 = ac(:,idx_C1)*100; results(4).C7p = sum(ac(:,idx_C7p),2)*100;

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
%  STYLE
%  ========================================================================
n_models = 4;
colors = {[0.00 0.45 0.74], ...   % blue   - Isothermal
          [0.85 0.33 0.10], ...   % red    - Haase
          [0.00 0.60 0.30], ...   % green  - Kempers
          [0.49 0.18 0.56]};      % purple - Firoozabadi

lstyles = {'-', '--', '-.', ':'};
lw      = 2.0;
lw_sat  = 1.4;
mk_sz   = 7;
GOC_field = 3647;
idx_plot = exp_h ~= h_ref;

ax_pos = [0.16 0.16 0.80 0.80];
y_lim  = [3580 3700];

fmt_ax = @(ax) set(ax, 'YDir', 'reverse', ...
    'FontName','Times New Roman', 'FontSize', 11, ...
    'TickDir','in', 'TickLength', [0.015 0.015], ...
    'XMinorTick','on', 'YMinorTick','on', 'LineWidth', 0.8);

%% ===== Figure 1: C1 mol% =====
fig1 = figure('Units','centimeters','Position',[2 2 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_lines = gobjects(n_models, 1);
for m = 1:n_models
    h_lines(m) = plot(results(m).C1, results(m).h, lstyles{m}, ...
        'Color', colors{m}, 'LineWidth', lw);
end
h_exp = plot(exp_C1(idx_plot), exp_h(idx_plot), 'ko', ...
    'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(ax); ylim(y_lim);
xlabel('C_{1} (mol%)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Depth (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

leg_str = cell(1, n_models);
for m = 1:n_models, leg_str{m} = results(m).name; end
lg = legend([h_exp; h_lines], [{'Measured (Pedersen 2006)'}, leg_str], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ===== Figure 2: C7+ mol% =====
fig2 = figure('Units','centimeters','Position',[4 3 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;
h_lines = gobjects(n_models, 1);
for m = 1:n_models
    h_lines(m) = plot(results(m).C7p, results(m).h, lstyles{m}, ...
        'Color', colors{m}, 'LineWidth', lw);
end
h_exp = plot(exp_C7p(idx_plot), exp_h(idx_plot), 'ko', ...
    'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(ax); ylim(y_lim);
xlabel('C_{7+} (mol%)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Depth (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

leg_str = cell(1, n_models);
for m = 1:n_models, leg_str{m} = results(m).name; end
lg = legend([h_exp; h_lines], [{'Measured (Pedersen 2006)'}, leg_str], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ===== Figure 3: Pressure & Psat (with GOC lines) =====
fig3 = figure('Units','centimeters','Position',[6 4 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_models, 1);
for m = 1:n_models
    h_lines(m) = plot(results(m).P/1e5, results(m).h, lstyles{m}, ...
        'Color', colors{m}, 'LineWidth', lw);
    v = ~isnan(results(m).Psat) & results(m).Psat > 0;
    plot(results(m).Psat(v)/1e5, results(m).h(v), lstyles{m}, ...
        'Color', colors{m}, 'LineWidth', lw_sat, 'HandleVisibility', 'off');
    yline(results(m).GOC, ':', 'Color', colors{m}, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
h_psat_dummy = plot(NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', lw_sat);
h_exp_P    = plot(exp_P(idx_plot), exp_h(idx_plot), 'ko', ...
    'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
h_exp_Psat = plot(exp_Psat(idx_plot), exp_h(idx_plot), 'k^', ...
    'MarkerSize', mk_sz, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
h_field_GOC = yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);

fmt_ax(ax); ylim(y_lim);
xlabel('Pressure (bar)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Depth (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

leg_str = cell(1, n_models);
for m = 1:n_models, leg_str{m} = sprintf('%s (GOC %.0f m)', results(m).name, results(m).GOC); end
lg = legend([h_exp_P; h_exp_Psat; h_field_GOC; h_psat_dummy; h_lines], ...
    [{'Measured P', 'Measured P_{sat}', sprintf('Field GOC (%d m)', GOC_field), ...
      'P_{sat} (model, thin)'}, leg_str], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% ===== Figure 4: Density =====
fig4 = figure('Units','centimeters','Position',[8 5 14 10],'Color','w','PaperPositionMode','auto');
ax = axes('Position', ax_pos); hold on; box on;

h_lines = gobjects(n_models, 1);
for m = 1:n_models
    v = ~isnan(results(m).rho) & results(m).rho > 0;
    h_lines(m) = plot(results(m).rho(v), results(m).h(v), lstyles{m}, ...
        'Color', colors{m}, 'LineWidth', lw);
end
h_exp = plot(exp_rho(idx_plot), exp_h(idx_plot), 'ko', ...
    'MarkerSize', mk_sz, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(ax); ylim(y_lim);
xlim(ax, [250 620]);
xlabel('Density (kg/m^{3})', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Depth (m)', 'FontSize', 12, 'FontName', 'Times New Roman');

leg_str = cell(1, n_models);
for m = 1:n_models, leg_str{m} = results(m).name; end
lg = legend([h_exp; h_lines], [{'Measured (Pedersen 2006)'}, leg_str], ...
    'Location', 'best', 'FontSize', 9, 'FontName', 'Times New Roman', ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
lg.ItemTokenSize = [22, 10];

%% Summary
fprintf('\n============ GOC Summary ============\n');
fprintf('%-15s  %8s\n', 'Model', 'GOC (m)');
fprintf('%s\n', repmat('-', 1, 26));
for m = 1:n_models
    fprintf('%-15s  %8.1f\n', results(m).name, results(m).GOC);
end
fprintf('%-15s  %8s\n', 'Experimental', '~3647');
fprintf('=====================================\n');