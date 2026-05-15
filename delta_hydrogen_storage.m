clear; clc; close all;

%% ========================================================================
%  DELTA SALT CAVERN — COMPOSITIONAL GRADING
%  Zhou et al. (2024), Table 1
%  30% H2, 120 bar, 320.43 K, 365.8 m thickness, h_ref = 1005.8 m
%  ========================================================================
h_ref     = 1005.8;
T_ref     = 320.43;
thickness = 365.8;
press_ref = 120e5;
R = 8.314462;
g = 9.80665;
dh = 1;
tol_sat = 1e-6;
maxiter_sat = 1500;
vt_method = 0;

dTdh = 0.023;

h_top = h_ref - thickness/2;
h_bot = h_ref + thickness/2;
h_down = (h_ref + dh : dh : h_bot)';
h_up   = (h_ref - dh : -dh : h_top)';

M_8   = [2.016; 28.014; 44.010; 16.043; 30.070; 44.096; 58.123; 58.123];
Tc_8  = [33.15; 126.20; 304.18; 190.60; 305.32; 369.83; 408.14; 425.12];
Pc_8  = [12.90; 33.98; 73.80; 46.00; 48.72; 42.48; 36.48; 37.96] * 1e5;
w_8   = [-0.22; 0.037; 0.239; 0.011; 0.0995; 0.1523; 0.1770; 0.2002];

BIP_8x8 = zeros(8);
BIP_8x8(1,2) = 0.10;  BIP_8x8(2,1) = 0.10;
BIP_8x8(1,3) = 0.10;  BIP_8x8(3,1) = 0.10;
BIP_8x8(1,4) = 0.04;  BIP_8x8(4,1) = 0.04;
BIP_8x8(1,5) = 0.08;  BIP_8x8(5,1) = 0.08;
BIP_8x8(1,6) = 0.10;  BIP_8x8(6,1) = 0.10;
BIP_8x8(1,7) = 0.10;  BIP_8x8(7,1) = 0.10;
BIP_8x8(1,8) = 0.11;  BIP_8x8(8,1) = 0.11;
BIP_8x8(2,3) = -0.017; BIP_8x8(3,2) = -0.017;
BIP_8x8(2,4) = 0.031;  BIP_8x8(4,2) = 0.031;
BIP_8x8(2,5) = 0.052;  BIP_8x8(5,2) = 0.052;
BIP_8x8(2,6) = 0.085;  BIP_8x8(6,2) = 0.085;
BIP_8x8(2,7) = 0.103;  BIP_8x8(7,2) = 0.103;
BIP_8x8(2,8) = 0.080;  BIP_8x8(8,2) = 0.080;
BIP_8x8(3,4) = 0.120;  BIP_8x8(4,3) = 0.120;
BIP_8x8(3,5) = 0.130;  BIP_8x8(5,3) = 0.130;
BIP_8x8(3,6) = 0.135;  BIP_8x8(6,3) = 0.135;
BIP_8x8(3,7) = 0.130;  BIP_8x8(7,3) = 0.130;
BIP_8x8(3,8) = 0.130;  BIP_8x8(8,3) = 0.130;
BIP_8x8(4,5) = 0.0026; BIP_8x8(5,4) = 0.0026;
BIP_8x8(4,6) = 0.014;  BIP_8x8(6,4) = 0.014;
BIP_8x8(4,7) = 0.016;  BIP_8x8(7,4) = 0.016;
BIP_8x8(4,8) = 0.020;  BIP_8x8(8,4) = 0.020;

Cp_8 = [
    27.14,   9.274e-3, -1.381e-5,  7.645e-9;   % H2
    31.15,  -0.014,     2.68e-5,  -1.17e-8;     % N2
    19.79,   0.073,    -5.60e-5,   1.72e-8;     % CO2
    19.25,   0.052,     1.20e-5,  -1.13e-8;     % CH4
     5.41,   0.178,    -6.94e-5,   8.71e-9;     % C2
    -4.22,   0.306,    -1.59e-4,   3.21e-8;     % C3
    -1.39,   0.385,    -1.85e-4,   2.90e-8;     % iC4
     9.49,   0.331,    -1.11e-4,  -2.82e-9;     % nC4
];

Hig_8 = [
    0;   % H2
    8008.5;                     % N2
    19023.5;                    % CO2
    0.0;                        % CH4
    10013.1;                    % C2
    20020.1;                    % C3
    28984.2;                    % iC4
    28984.2;                    % nC4
];

cases = struct();

cases(1).name = 'Pure CH_4';
cases(1).comp = [0.30; 0.70];
cases(1).M    = M_8([1 4]);
cases(1).Tc   = Tc_8([1 4]);
cases(1).Pc   = Pc_8([1 4]);
cases(1).w    = w_8([1 4]);
cases(1).BIP  = [0 0.04; 0.04 0];
cases(1).Cp   = Cp_8([1 4], :);
cases(1).Hig  = Hig_8([1 4]);

cases(2).name = 'Rich NG';
cases(2).comp = [0.30; 0.01; 0.015; 0.25; 0.305; 0.064; 0.029; 0.027];
cases(2).M = M_8; cases(2).Tc = Tc_8; cases(2).Pc = Pc_8; cases(2).w = w_8;
cases(2).BIP = BIP_8x8; cases(2).Cp = Cp_8; cases(2).Hig = Hig_8;

cases(3).name = 'Pure C_2H_6';
cases(3).comp = [0.30; 0.70];
cases(3).M    = M_8([1 5]);
cases(3).Tc   = Tc_8([1 5]);
cases(3).Pc   = Pc_8([1 5]);
cases(3).w    = w_8([1 5]);
cases(3).BIP  = [0 0.08; 0.08 0];
cases(3).Cp   = Cp_8([1 5], :);
cases(3).Hig  = Hig_8([1 5]);

cases(4).name = 'Pure CO_2';
cases(4).comp = [0.30; 0.70];
cases(4).M    = M_8([1 3]);
cases(4).Tc   = Tc_8([1 3]);
cases(4).Pc   = Pc_8([1 3]);
cases(4).w    = w_8([1 3]);
cases(4).BIP  = [0 0.10; 0.10 0];
cases(4).Cp   = Cp_8([1 3], :);
cases(4).Hig  = Hig_8([1 3]);

ncases = length(cases);
colors = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.30 0.75 0.93], [0.49 0.18 0.56]};
lstyles = {'-', '--', '-.', ':'};

%% ========================================================================
%  PART 1: COMPOSITIONAL GRADING
%  ========================================================================
if dTdh == 0
    fprintf('=== DELTA: CUSHION GAS COMPARISON (isothermal, 120 bar, %.2f K) ===\n', T_ref);
else
    fprintf('=== DELTA: CUSHION GAS COMPARISON (dT/dh = %.3f K/m, 120 bar, %.2f K) ===\n', dTdh, T_ref);
end

results = struct();
for c = 1:ncases
    comp_c = cases(c).comp;
    idx = comp_c == 0;
    comp_c(idx) = []; M_c = cases(c).M; M_c(idx) = [];
    Tc_c = cases(c).Tc; Tc_c(idx) = []; Pc_c = cases(c).Pc; Pc_c(idx) = [];
    w_c = cases(c).w; w_c(idx) = []; BIP_c = cases(c).BIP; BIP_c(idx,:) = []; BIP_c(:,idx) = [];
    Cp_c = cases(c).Cp; Cp_c(idx,:) = []; Hig_c = cases(c).Hig; Hig_c(idx) = [];
    comp_c = comp_c / sum(comp_c);

    [all_h, all_z, all_P] = run_hasse(comp_c, press_ref, T_ref, ...
        Pc_c, Tc_c, w_c, M_c, BIP_c, Cp_c, Hig_c, h_ref, h_down, h_up, dTdh, vt_method);

    results(c).h = all_h;
    results(c).z = all_z;
    results(c).P = all_P;
    results(c).xH2 = all_z(:,1)*100;
    results(c).name = cases(c).name;
    fprintf('  %15s: x_H2 = %.4f to %.4f mol%%  (Delta = %.4f)\n', ...
        cases(c).name, results(c).xH2(1), results(c).xH2(end), results(c).xH2(1)-results(c).xH2(end));
end

figure('Position', [100 100 750 600]);
hold on;
for c = 1:ncases
    plot(results(c).xH2, results(c).h, 'Color', colors{c}, 'LineStyle', lstyles{c}, 'LineWidth', 1.8);
end
xlabel('H_2 (mol%)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse', 'FontSize', 11);
legend({results.name}, 'Location', 'best', 'FontSize', 10);
grid on;
if dTdh == 0
    title('Delta Salt Cavern — Effect of Cushion Gas Type (isothermal, 120 bar)', 'FontSize', 12);
else
    title(sprintf('Delta Salt Cavern — Effect of Cushion Gas Type (dT/dh = %.3f K/m, 120 bar)', dTdh), 'FontSize', 12);
end

%% ========================================================================
%  PART 2: SATURATION PRESSURE vs RESERVOIR PRESSURE — Rich NG
%  ========================================================================
fprintf('\n=== RICH NG: SATURATION PRESSURE vs DEPTH ===\n');

c = 2;
comp_c = cases(c).comp;
idx = comp_c == 0;
comp_c(idx) = []; M_c = cases(c).M; M_c(idx) = [];
Tc_c = cases(c).Tc; Tc_c(idx) = []; Pc_c = cases(c).Pc; Pc_c(idx) = [];
w_c = cases(c).w; w_c(idx) = []; BIP_c = cases(c).BIP; BIP_c(idx,:) = []; BIP_c(:,idx) = [];
nc_c = length(comp_c);
components_c = cell(nc_c, 1);

sample_step = 10;
all_h_c = results(c).h;
idx_sample = 1 : sample_step : length(all_h_c);
if idx_sample(end) ~= length(all_h_c)
    idx_sample = [idx_sample, length(all_h_c)];
end

h_sample = all_h_c(idx_sample);
P_sample = results(c).P(idx_sample);
z_sample = results(c).z(idx_sample, :);

Pdew_sample = nan(length(idx_sample), 1);
for k = 1:length(idx_sample)
    z_k = z_sample(k,:)';
    z_k = max(z_k, 1e-15);
    z_k = z_k / sum(z_k);
    P_k = P_sample(k);
    try
        [Pdew, ~] = pressdew_multicomp_newton(z_k, P_k, T_ref, Pc_c, Tc_c, w_c, BIP_c, tol_sat, maxiter_sat, 'PR', components_c);
        if isreal(Pdew) && Pdew > 0 && isfinite(Pdew)
            Pdew_sample(k) = Pdew;
        end
    catch
    end
    if isnan(Pdew_sample(k))
        try
            [Pdew, ~] = pressdew_multicomp_ss(z_k, P_k, T_ref, Pc_c, Tc_c, w_c, BIP_c, tol_sat, maxiter_sat, 'PR', components_c);
            if isreal(Pdew) && Pdew > 0 && isfinite(Pdew)
                Pdew_sample(k) = Pdew;
            end
        catch
        end
    end
end

figure('Position', [100 100 900 650]);
hold on;
plot_handles = [];
legend_entries = {};

h1 = plot(P_sample/1e5, h_sample, 'Color', [0.47 0.67 0.19], 'LineStyle', '-', 'LineWidth', 1.8);
plot_handles(end+1) = h1;
legend_entries{end+1} = 'P_{res}';

valid = ~isnan(Pdew_sample);
if any(valid)
    h2 = plot(Pdew_sample(valid)/1e5, h_sample(valid), 'Color', [0.64 0.08 0.18], 'LineStyle', '--', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
    plot_handles(end+1) = h2;
    legend_entries{end+1} = 'P_{dew}';
    margin = min((P_sample(valid) - Pdew_sample(valid)) ./ P_sample(valid) * 100);
else
    margin = Inf;
end
fprintf('  Rich NG: min margin to dew point = %.1f%%\n', margin);

xlabel('Pressure (bar)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse', 'FontSize', 11);
legend(plot_handles, legend_entries, 'Location', 'best', 'FontSize', 10);
grid on;
title('Rich NG — Reservoir Pressure vs Dew Point (120 bar, 320.43 K)', 'FontSize', 12);

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================
function [all_h, all_z, all_P] = run_hasse(comp_ref, P_ref, T_ref, Pc, Tc, w, M_gmol, BIP, Cp_coeffs, H_ig_ref, h_ref, h_down, h_up, dTdh, vt_method)

    nc = length(comp_ref);
    n_down = length(h_down);
    n_up   = length(h_up);

    comp_down = zeros(n_down, nc); P_down = zeros(n_down, 1);
    comp_up   = zeros(n_up, nc);   P_up   = zeros(n_up, 1);

    h_prev = h_ref; comp_prev = comp_ref(:); P_prev = P_ref; T_prev = T_ref;
    for k = 1:n_down
        [comp_next, P_next, ~, ~, ~] = main_hasse(h_down(k), h_prev, comp_prev, P_prev, T_prev, dTdh, ...
            Pc, Tc, w, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method);
        comp_down(k,:) = comp_next(:)';
        P_down(k) = P_next;
        T_next = T_prev + dTdh * (h_down(k) - h_prev);
        h_prev = h_down(k); comp_prev = comp_next; P_prev = P_next; T_prev = T_next;
    end

    h_prev = h_ref; comp_prev = comp_ref(:); P_prev = P_ref; T_prev = T_ref;
    for k = 1:n_up
        [comp_next, P_next, ~, ~, ~] = main_hasse(h_up(k), h_prev, comp_prev, P_prev, T_prev, dTdh, ...
            Pc, Tc, w, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method);
        comp_up(k,:) = comp_next(:)';
        P_up(k) = P_next;
        T_next = T_prev + dTdh * (h_up(k) - h_prev);
        h_prev = h_up(k); comp_prev = comp_next; P_prev = P_next; T_prev = T_next;
    end

    all_h = [flipud(h_up); h_ref; h_down];
    all_z = [flipud(comp_up); comp_ref(:)'; comp_down];
    all_P = [flipud(P_up); P_ref; P_down];
end