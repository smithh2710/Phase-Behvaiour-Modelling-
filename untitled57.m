clear; clc; close all;

%  DELTA SALT CAVERN — H2/CH4 BINARY: TEMPERATURE GRADIENT SENSITIVITY
%  30% H2 / 70% CH4, 120 bar, 320.43 K, 365.8 m thickness

h_ref     = 1005.8;
T_ref     = 320.43;
thickness = 365.8;
press_ref = 120e5;
R = 8.314462;
g = 9.80665;
dh = 5;

h_top = h_ref - thickness/2;
h_bot = h_ref + thickness/2;
h_down = (h_ref + dh : dh : h_bot)';
h_up   = (h_ref - dh : -dh : h_top)';

comp = [0.30; 0.70];
M     = [2.016; 16.043];
Tc    = [33.20; 190.60];
Pc    = [12.97; 46.00] * 1e5;
w     = [-0.22; 0.0080];
BIP   = [0 0.0156; 0.0156 0];

Cp_coeffs = [
    27.14,   9.274e-3, -1.381e-5,  7.645e-9;    % H2  (RPP)
    19.25,   0.052,     1.20e-5,  -1.13e-8;     % CH4 (Pedersen)
];

H_ig_ref = [
    0.0;   % H2  — Pedersen formula extrapolation
    0.0;   % CH4 (Pedersen 2006 Table 3)
];

nc = length(comp);
vt_method = 0;

%% ========================================================================
%  TEMPERATURE GRADIENT SENSITIVITY
%  dTdh in K/m — positive means T increases with depth (normal geothermal)
%  Typical salt cavern: 0.01–0.02 K/m (salt has high conductivity)
%  Typical reservoir:   0.025–0.035 K/m
%  ========================================================================

dTdh_vals = [0, 0.01, 0.023, 0.03, 0.05];
n_cases = length(dTdh_vals);

colors = {[0.00 0.00 0.00], [0.00 0.45 0.74], [0.85 0.33 0.10], ...
          [0.93 0.69 0.13], [0.47 0.67 0.19], [0.64 0.08 0.18]};
lstyles = {'-', '-', '--', '-.', ':', '-'};
lwidths = [2.5, 1.5, 1.5, 1.5, 1.5, 1.5];

fprintf('=== H2/CH4 BINARY: TEMPERATURE GRADIENT SENSITIVITY (120 bar, 320.43 K) ===\n');
fprintf('%-12s %-12s %-12s %-12s %-14s\n', 'dT/dh', 'x_H2 top', 'x_H2 bot', 'Delta', 'T range (K)');
fprintf('%s\n', repmat('-', 1, 65));

results = struct();
for c = 1:n_cases
    dTdh = dTdh_vals(c);

    [all_h, all_z, all_P] = run_hasse(comp, press_ref, T_ref, ...
        Pc, Tc, w, M, BIP, Cp_coeffs, H_ig_ref, h_ref, h_down, h_up, dTdh, vt_method);

    results(c).h    = all_h;
    results(c).z    = all_z;
    results(c).P    = all_P;
    results(c).xH2  = all_z(:,1) * 100;
    results(c).dTdh = dTdh;

    T_top = T_ref + dTdh * (h_top - h_ref);
    T_bot = T_ref + dTdh * (h_bot - h_ref);

    if dTdh == 0
        results(c).name = 'Isothermal';
    else
        results(c).name = sprintf('dT/dh = %.3f K/m', dTdh);
    end

    fprintf('%-12.3f %-12.4f %-12.4f %-12.4f %-7.1f – %-7.1f\n', ...
        dTdh, results(c).xH2(1), results(c).xH2(end), ...
        results(c).xH2(1) - results(c).xH2(end), T_top, T_bot);
end

%% ========================================================================
%  FIGURE 1: H2 mol% vs Depth — all dT/dh cases
%  ========================================================================
figure('Position', [100 100 800 650]);
hold on;
for c = 1:n_cases
    plot(results(c).xH2, results(c).h, ...
        'Color', colors{c}, 'LineStyle', lstyles{c}, 'LineWidth', lwidths(c));
end
xlabel('H_2 (mol%)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse', 'FontSize', 11);
legend({results.name}, 'Location', 'best', 'FontSize', 9);
grid on;
title('H_2/CH_4 Binary — Effect of Temperature Gradient on H_2 Segregation (120 bar)', 'FontSize', 12);

%% ========================================================================
%  FIGURE 2: Delta x_H2 (top-bottom) vs dT/dh
%  ========================================================================
delta_xH2 = zeros(n_cases, 1);
for c = 1:n_cases
    delta_xH2(c) = results(c).xH2(1) - results(c).xH2(end);
end

figure('Position', [100 100 600 400]);
plot(dTdh_vals, delta_xH2, 'o-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.00 0.45 0.74]);
xlabel('dT/dh (K/m)'); ylabel('\Delta x_{H_2} (mol%)');
set(gca, 'FontSize', 11);
grid on;
title('H_2 Segregation Magnitude vs Temperature Gradient', 'FontSize', 12);

%% ========================================================================
%  LOCAL: sequential march using main_hasse
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