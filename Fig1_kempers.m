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
press_ref = 380.20e5;
temp_ref = 138.70 + 273.15;
dTdh = 0.0260;
R = 8.3144598;

Zc = 0.29056 - 0.08775 .* acentric;
for i = 2:length(Zc)
    if Zc(i) >= Zc(i-1)
        Zc(i) = Zc(i-1) - 0.008;
    end
end
Vc = Zc .* R .* Tc ./ Pc;

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

tol = 1e-10;
maxiter = 500;

vt_method = 7;
vt_params.c_custom = zeros(n,1);
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

%% GOC detection
[GOC_h, ~, ~, ~] = detect_sgoc(2, comp_ref, press_ref, temp_ref, h_ref, dTdh, ...
    [3580 3682], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
fprintf('GOC: %.1f m\n', GOC_h);
GOC_T = temp_ref + dTdh * (GOC_h - h_ref);

%% Build depth arrays relative to actual GOC
measured_depths = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];

oil_measured = measured_depths(measured_depths > GOC_h + 2);
h_oil = unique([GOC_h; GOC_h+5; GOC_h+10; oil_measured]);
h_oil = sort(h_oil);

gas_measured = measured_depths(measured_depths < GOC_h - 2);
h_gas = unique([GOC_h; GOC_h-5; GOC_h-10; gas_measured]);
h_gas = sort(h_gas, 'descend');

n_oil = length(h_oil);
n_gas = length(h_gas);

%% Oil leg
oil_comp = zeros(n_oil, n);
oil_P    = NaN(n_oil,1); oil_T    = NaN(n_oil,1);
oil_Pbub = NaN(n_oil,1); oil_rho  = NaN(n_oil,1);

for i = 1:n_oil
    try
        [ch, Ph, Th, Pb, ~] = main_kempers(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        oil_comp(i,:) = ch(:)';
        oil_P(i) = Ph; oil_T(i) = Th; oil_Pbub(i) = Pb;
        oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch ME
        fprintf('Oil error at h=%.1f: %s\n', h_oil(i), ME.message);
    end
end

%% Flash at GOC to get equilibrium vapor
[~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), oil_T(1), ...
    Pc, Tc, acentric, BIP, tol, maxiter);

%% Gas leg
gas_comp = zeros(n_gas, n);
gas_P    = NaN(n_gas,1); gas_T    = NaN(n_gas,1);
gas_Pdew = NaN(n_gas,1); gas_rho  = NaN(n_gas,1);

for i = 1:n_gas
    try
        [ch, Ph, Th, ~, Pd] = main_kempers(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), oil_T(1), dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        gas_comp(i,:) = ch(:)';
        gas_P(i) = Ph; gas_T(i) = Th; gas_Pdew(i) = Pd;
        if isnan(Pd) && all(ch > 0)
            try
                [Pd_ss, ~] = pressdew_multicomp_ss(ch(:), Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
                gas_Pdew(i) = Pd_ss;
            catch
            end
        end
        gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch ME
        fprintf('Gas error at h=%.1f: %s\n', h_gas(i), ME.message);
    end
end

gas_Pdew(1) = oil_Pbub(1);

oil_C1  = oil_comp(:,3) * 100;
gas_C1  = gas_comp(:,3) * 100;
oil_C7p = sum(oil_comp(:,11:end), 2) * 100;
gas_C7p = sum(gas_comp(:,11:end), 2) * 100;

%% Experimental data (Pedersen & Hjermstad 2006)
exp_h       = [3638.2; 3644.3; 3651.1; 3661.6; 3676.0; 3682.8];
exp_Pres    = [377.8; 377.9; 378.2; 378.8; 379.6; 380.2];
exp_Psat    = [375.5; 372.8; 364.2; 364.5; 360.1; 353.8];
exp_rho     = [367;    376;    574;    581;    595;    601];
exp_C1      = [68.861; 68.546; 56.142; 55.261; 54.253; 53.871];
exp_c7p     = [1.482+1.595+1.031+5.231; 1.519+1.610+1.048+5.080; ...
               2.474+2.583+1.695+13.78; 2.520+2.667+1.779+14.491; ...
               2.579+2.777+1.869+15.282; 2.630+2.823+1.897+15.783];

%% Plot
fn  = 'Times New Roman';
fsz = 10;
lw  = 1.8;
col = [0.00 0.45 0.74];
mk  = [0.6 0.6 0.6];

figure('Units','centimeters','Position',[2 2 28 20],'Color','w','PaperPositionMode','auto');

% (a) Pressure
subplot(2,2,1); hold on;
all_h_p    = [flipud(h_gas); h_oil];
all_Pres   = [flipud(gas_P); oil_P];
all_Psat   = [flipud(gas_Pdew); oil_Pbub];
h1 = plot(all_Pres/1e5, all_h_p, '-', 'Color', col, 'LineWidth', lw);
h2 = plot(all_Psat/1e5, all_h_p, '--', 'Color', col, 'LineWidth', lw);
h3 = plot(exp_Pres, exp_h, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', mk, 'LineWidth', 0.8);
h4 = plot(exp_Psat, exp_h, 'k^', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
set(gca,'YDir','reverse','FontSize',fsz,'FontName',fn,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','LineWidth',0.6,'Box','on');
ylim([3620 3700]);
xlabel('Pressure (bar)','FontSize',fsz,'FontName',fn);
ylabel('Depth (m)','FontSize',fsz,'FontName',fn);
title('(a)','FontSize',fsz+1,'FontWeight','bold','FontName',fn);
legend([h1 h2 h3 h4], {'P_{res}','P_{sat}','Measured P_{res}','Measured P_{sat}'}, ...
    'Location','best','FontSize',8,'FontName',fn,'Box','off');

% (b) Density
subplot(2,2,2); hold on;
h1 = plot(oil_rho, h_oil, '-', 'Color', col, 'LineWidth', lw);
plot(gas_rho, h_gas, '-', 'Color', col, 'LineWidth', lw);
if ~isnan(gas_rho(1)) && ~isnan(oil_rho(1))
    plot([gas_rho(1), oil_rho(1)], [GOC_h, GOC_h], '-', 'Color', col, 'LineWidth', lw);
end
h2 = plot(exp_rho, exp_h, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', mk, 'LineWidth', 0.8);
set(gca,'YDir','reverse','FontSize',fsz,'FontName',fn,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','LineWidth',0.6,'Box','on');
xlim([300 650]); ylim([3620 3700]);
xlabel('Density (kg/m^3)','FontSize',fsz,'FontName',fn);
ylabel('Depth (m)','FontSize',fsz,'FontName',fn);
title('(b)','FontSize',fsz+1,'FontWeight','bold','FontName',fn);
legend([h1 h2], {'Kempers, no VT','Measured (Pedersen 2006)'}, ...
    'Location','southeast','FontSize',8,'FontName',fn,'Box','off');

% (c) C1
subplot(2,2,3); hold on;
v_o = oil_C1 > 0; v_g = gas_C1 > 0;
h1 = plot(oil_C1(v_o), h_oil(v_o), '-', 'Color', col, 'LineWidth', lw);
plot(gas_C1(v_g), h_gas(v_g), '-', 'Color', col, 'LineWidth', lw);
if gas_C1(1) > 0 && oil_C1(1) > 0
    plot([gas_C1(1), oil_C1(1)], [GOC_h, GOC_h], '-', 'Color', col, 'LineWidth', lw);
end
h2 = plot(exp_C1, exp_h, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', mk, 'LineWidth', 0.8);
set(gca,'YDir','reverse','FontSize',fsz,'FontName',fn,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','LineWidth',0.6,'Box','on');
xlim([50 72]); ylim([3620 3700]);
xlabel('C_1 (mol%)','FontSize',fsz,'FontName',fn);
ylabel('Depth (m)','FontSize',fsz,'FontName',fn);
title('(c)','FontSize',fsz+1,'FontWeight','bold','FontName',fn);
legend([h1 h2], {'Kempers, no VT','Measured (Pedersen 2006)'}, ...
    'Location','southeast','FontSize',8,'FontName',fn,'Box','off');

% (d) C7+
subplot(2,2,4); hold on;
v_o = oil_C7p > 0; v_g = gas_C7p > 0;
h1 = plot(oil_C7p(v_o), h_oil(v_o), '-', 'Color', col, 'LineWidth', lw);
plot(gas_C7p(v_g), h_gas(v_g), '-', 'Color', col, 'LineWidth', lw);
if gas_C7p(1) > 0 && oil_C7p(1) > 0
    plot([gas_C7p(1), oil_C7p(1)], [GOC_h, GOC_h], '-', 'Color', col, 'LineWidth', lw);
end
h2 = plot(exp_c7p, exp_h, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', mk, 'LineWidth', 0.8);
set(gca,'YDir','reverse','FontSize',fsz,'FontName',fn,'TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','LineWidth',0.6,'Box','on');
ylim([3620 3700]);
xlabel('C_{7+} (mol%)','FontSize',fsz,'FontName',fn);
ylabel('Depth (m)','FontSize',fsz,'FontName',fn);
title('(d)','FontSize',fsz+1,'FontWeight','bold','FontName',fn);
legend([h1 h2], {'Kempers, no VT','Measured (Pedersen 2006)'}, ...
    'Location','southeast','FontSize',8,'FontName',fn,'Box','off');