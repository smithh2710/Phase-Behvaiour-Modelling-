
clc ; clear all ; close all ; 

components = {'N2', 'CO2', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', 'C6', 'C7', 'C8', 'C9', 'C10-C11', 'C12-C13', 'C14-C16','C17-C18', 'C19-C21', 'C22-C24', 'C25-C29', 'C30-C37', 'C38-C80'};
comp_ref = [0.42; 0.71; 49.885; 7.771; 6.741; 1.030; 3.21; 1.16; 1.55; 1.88 ;3.490; 3.730 ; 2.31 ; 3.293; 2.620; 2.961; 1.480; 1.673; 1.187; 1.264;0.982; 0.652];
comp_ref = comp_ref ./ sum(comp_ref) ;
M_gmol = [28.01,44.01 ,16.04 ,30.07 ,44.10 ,58.12 ,58.12 ,72.15 ,72.15 ,86.18 ,96 ,107 ,121 ,140.09 ,167.57 ,204.75 ,243.6 ,275.27 ,317.02 ,370.39 ,456.83 ,640.76]';
acentric = [0.0400 ,0.2250 ,0.0080 ,0.0980 ,0.1520 ,0.1760 ,0.1930 ,0.2270 ,0.2510 ,0.2960 ,0.3375 ,0.3743 ,0.4204 ,0.4833 ,0.5703 ,0.6847 ,0.7947 ,0.8805 ,0.9836 ,1.0983 ,1.2292 ,1.1596]' ;
Tc = [-146.95 ,31.05 ,-82.55 ,32.25 ,96.65 ,134.95 ,152.05 ,187.25 ,196.45 234.25 ,301.63 ,323.97 ,349.62 ,381.84 ,422.96 ,473.01 ,519.43 ,555.64 ,600.26 ,654.86 ,737.85 ,912.38]' + 273.15 ;
Pc = [33.94 ,73.76 ,46.00 ,48.84 ,42.46 ,36.48 ,38 ,33.84 ,33.74 ,29.69 ,29.59 ,27.37 ,25.01 ,22.59 ,20.12 ,17.91 ,16.4 ,15.55 ,14.71 ,13.96 ,13.16 ,12.27 ]'  * 1e5 ; 
c_JY = [-4.23; -1.64; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12;1.39; 15.53 ; 19.96 ; 25.14 ; 31.89 ; 39.65 ; 46.52 ; 49.62 ; 50.37 ; 48.37; 42.77 ; 27.92; -6.04];

h_ref = 204;
press_ref = 286e5;
temp_ref = 94 + 273.15;

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
for i = 10:n
    BIP(1,i) = 0.0800; BIP(i,1) = 0.0800;
end
for i = 3:10
    BIP(2,i) = 0.1200; BIP(i,2) = 0.1200;
end
for i = 11:n
    BIP(2,i) = 0.100; BIP(i,2) = 0.100;
end

Vc = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0;432; 492; 548; 600.22; 697.08; 810.90; 910.91; 980.42; 1061.51;1152.43; 1287.26; 1648.29];
Zc = Pc .* (Vc * 1e-6) ./ (R .* Tc);

dTdh = 0.025;
tol = 1e-10;
maxiter = 500;

Cp_coeffs = [
%   c1          c2        c3          c4
  31.15,    -0.014,    2.68e-5,   -1.17e-8;   % N2
  19.79,     0.073,   -5.60e-5,    1.72e-8;   % CO2
  19.25,     0.052,    1.20e-5,   -1.13e-8;   % C1
   5.41,     0.178,   -6.94e-5,    8.71e-9;   % C2
  -4.22,     0.306,   -1.59e-4,    3.21e-8;   % C3
  -1.39,     0.385,   -1.85e-4,    2.90e-8;   % iC4
   9.49,     0.331,   -1.11e-4,   -2.82e-9;   % nC4
  -9.52,     0.507,   -2.73e-4,    5.73e-8;   % iC5
  -3.63,     0.487,   -2.73e-4,    5.30e-8;   % nC5
  -4.41,     0.582,   -3.12e-4,    6.49e-8;   % C6
  -5.15,     0.676,   -3.65e-4,    7.66e-8;   % nC7
  -6.30,     0.771,   -4.20e-4,    8.85e-8;   % nC8
   3.14,     0.677,   -1.93e-4,   -2.98e-8;   % nC9
  25.20,     0.830,   -3.23e-4,    4.05e-8;   % C10-C11
  30.14,     0.993,   -3.87e-4,    4.86e-8;   % C12-C13
  36.83,     1.213,   -4.72e-4,    5.93e-8;   % C14-C16
  43.82,     1.443,   -5.62e-4,    7.06e-8;   % C17-C18
  49.52,     1.630,   -6.35e-4,    7.98e-8;   % C19-C21
  57.03,     1.878,   -7.31e-4,    9.19e-8;   % C22-C24
  66.63,     2.194,   -8.55e-4,    1.07e-7;   % C25-C29
  82.18,     2.706,   -1.05e-3,    1.32e-7;   % C30-C37
 115.26,     3.795,   -1.48e-3,    1.86e-7;   % C38-C80
];

% H_ig_ref = [8008.5; 19023.5; 0.0; 10013.1; 20020.1; 28984.2; 28984.2; 39079.8;29054.7; 49022.7; 55626.4; 63278.8; 75018.2; 86305.5; 105436.4;131245.8; 158307.4; 180360.2; 209369.7; 246518.5; 306624.4; 434697.4];

H_ig_ref = [8330.789; 19459.101; 2.642; 9761.134; 19519.622; 29278.121; 29278.121; 39036.609; 39036.609; 48795.1026; ...
            58553.595; 68312.084; 78069.882; 86301.327; 105418.986; 131284.869; 158312.556; 180345.162; 209390.365; ...
            246519.548; 306655.264; 434614.185];


% H_ig_ref_per_mass = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33; 2; 2; 2; 2; 93; 93; 93; 31; 31; 31; 8; 8];
% H_ig_ref = H_ig_ref_per_mass  .* M_gmol * R ;

vt_params.c_custom = c_JY;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;
vt_method = 5 ;

 [GOC_depth, GOC_pressure, ~, GOC_temp] = detect_sgoc(1, comp_ref, press_ref, temp_ref, h_ref, dTdh, [20 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)
[comp_h, P_h, T_h, Pb,~] = main_hasse(327, h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)

comp_h = comp_h*100 ;  
inorg = sum(comp_h(1:2))
c1 = comp_h(3)
c2_c6 = sum(comp_h(4:10) ) 
c7_plus = sum(comp_h(11:n))

[GOR, Bo] = calculate_GOR_STO(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol,'vt_method', vt_method, 'vt_params', vt_params)
[rho, ~, ~, ~] = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params)
%% Oil zone
h_oil = [GOC_depth,GOC_depth + 10 , 175, 204, 228, 327];
n_oil = length(h_oil);
oil_comp = zeros(n_oil, n);
oil_P    = zeros(n_oil, 1);
oil_T    = zeros(n_oil, 1);
oil_Pbub = zeros(n_oil, 1);
oil_rho  = zeros(n_oil, 1);
oil_GOR  = zeros(n_oil, 1);

for i = 1:n_oil
    [comp_h, P_h, T_h, Pb, ~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    oil_comp(i,:) = comp_h';
    oil_P(i)    = P_h;
    oil_T(i)    = T_h;
    oil_Pbub(i) = Pb;
    try
        oil_rho(i) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch
        oil_rho(i) = NaN;
    end
    try
        oil_GOR(i) = calculate_GOR_STO(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, 'vt_method', vt_method, 'vt_params', vt_params);
    catch
        oil_GOR(i) = NaN;
    end
end

%% GOC -> vapor composition
[~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), oil_T(1), Pc, Tc, acentric, BIP, tol, maxiter);

%% Gas zone
h_gas = [GOC_depth,(GOC_depth-30), 100, 80, 60, 40, 20, 0];
n_gas = length(h_gas);
gas_comp = zeros(n_gas, n);
gas_P    = zeros(n_gas, 1);
gas_T    = zeros(n_gas, 1);
gas_Pdew = zeros(n_gas, 1);
gas_rho  = zeros(n_gas, 1);
gas_GOR  = zeros(n_gas, 1);

for i = 1:n_gas
    [comp_h, P_h, T_h, ~, Pd] = main_hasse(h_gas(i), GOC_depth, comp_gas_ref, oil_P(1), oil_T(1), dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    gas_comp(i,:) = comp_h';
    gas_P(i)    = P_h;
    gas_T(i)    = T_h;
    gas_Pdew(i) = Pd;
    try
        gas_rho(i) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch
        gas_rho(i) = NaN;
    end
    try
        gas_GOR(i) = calculate_GOR_STO(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, 'vt_method', vt_method, 'vt_params', vt_params);
    catch
        gas_GOR(i) = NaN;
    end
end

%% Combine (gas top -> GOC_gas, GOC_oil -> bottom)
all_h    = [flipud(h_gas'); h_oil'];
all_P    = [flipud(gas_P);  oil_P];
all_Psat = [flipud(gas_Pdew); oil_Pbub];
all_rho  = [flipud(gas_rho);  oil_rho];
all_GOR  = [flipud(gas_GOR);  oil_GOR];


exp_h    = [0; 175; 204; 228; 327];
exp_P    = [279; 284; 286; 287; 293];
exp_Psat = [270; 272; 267; 265; 242];

figure('Position', [50 100 500 450], 'Color', 'w'); hold on;
plot(all_P/1e5, all_h, 'b-', 'LineWidth', 2);
plot(all_Psat/1e5, all_h, 'b--', 'LineWidth', 2);
plot(exp_P, exp_h, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(exp_Psat, exp_h, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
yline(GOC_depth, 'r:', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('Pressure [bar]'); ylabel('Depth [m]');
title('Reservoir 1 — Haase Model, No VT');
legend('P_{res}', 'P_{sat}', 'Exp P_{res}', 'Exp P_{sat}', 'GOC', 'Location', 'southeast');
grid on;
 
% Plot 2: Density vs Depth
figure('Position', [600 100 500 450], 'Color', 'w'); hold on;
plot(all_rho, all_h, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
yline(GOC_depth, 'r:', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('Density [kg/m³]'); ylabel('Depth [m]');
title('Reservoir 1 — Haase Model, Density vs Depth');
legend('Density', 'GOC', 'Location', 'best');
grid on;

% Plot 3: GOR vs Depth
figure('Position', [1150 100 500 450], 'Color', 'w'); hold on;
plot(all_GOR, all_h, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
yline(GOC_depth, 'r:', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('GOR [Sm³/Sm³]'); ylabel('Depth [m]');
title('Reservoir 1 — Haase Model, GOR vs Depth');
legend('GOR', 'GOC', 'Location', 'best');
grid on;