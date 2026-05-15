clear; clc; close all;


%  DELTA SALT CAVERN — HEAVY NG CUSHION GAS
%  30% H2 + Heavy NG, 120 bar, 320.43 K, 365.8 m thickness

Components =[ 'H2','N2', 'CO2', 'CH4', 'C2H6', 'C3H8', 'iC4H10', 'nC4H10' ] ;
h_ref     = 1005.8;
temp      = 320.43;
thickness = 365.8;
R = 8.314462;
g = 9.80665;
dh = 1;
tol = 1e-12;
maxiter = 300;
press_ref = 180e5;

h_top = h_ref - thickness/2;
h_bot = h_ref + thickness/2;
h_down = (h_ref + dh : dh : h_bot)';
h_up   = (h_ref - dh : -dh : h_top)';

comp = [0.30; 0.01; 0.015; 0.25; 0.18; 0.13; 0.06; 0.055];
M     = [2.016; 28.014; 44.010; 16.043; 30.070; 44.096; 58.124; 58.124];
Tc    = [33.20; 126.20; 304.2; 190.60; 305.400; 369.80; 408.1; 425.2];
Pc    = [12.97; 33.94; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00] * 1e5;
w     = [-0.22; 0.04; 0.2250; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930];

BIP = zeros(8);
BIP(1,2) = 0.0738;  BIP(2,1) = 0.0738;
BIP(1,3) = -0.1622;  BIP(3,1) = -0.1622;
BIP(1,4) = 0.0156;  BIP(4,1) = 0.0156;
BIP(1,5) = -0.0667;  BIP(5,1) = -0.0667;
BIP(1,6) = -0.0833;  BIP(6,1) = -0.0833;
BIP(1,7) = 0.1127;  BIP(7,1) = 0.1127;
BIP(1,8) = -0.3970;  BIP(8,1) = -0.3970;

BIP(2,3) = -0.017; BIP(3,2) = -0.017;
BIP(2,4) = 0.031;  BIP(4,2) = 0.031;
BIP(2,5) = 0.0515;  BIP(5,2) = 0.0515;
BIP(2,6) = 0.0852;  BIP(6,2) = 0.0852;
BIP(2,7) = 0.1033;  BIP(7,2) = 0.1033;
BIP(2,8) = 0.080;  BIP(8,2) = 0.080;

BIP(3,4) = 0.120;  BIP(4,3) = 0.120;
BIP(3,5) = 0.120;  BIP(5,3) = 0.120;
BIP(3,6) = 0.120;  BIP(6,3) = 0.120;
BIP(3,7) = 0.120;  BIP(7,3) = 0.120;
BIP(3,8) = 0.120;  BIP(8,3) = 0.120;

nc = length(comp);
Cp_dummy = zeros(nc, 5);
Hig_dummy = zeros(nc, 1);
vt_method = 0;
vt_params = struct('components', {cell(nc, 1)});
depth_range = [h_top, h_bot];


%% ========================================================================
%  PART 1: COMPOSITIONAL GRADING
%  ========================================================================
fprintf('=== HEAVY NG: COMPOSITIONAL GRADING (180 bar, 320.43 K) ===\n');

[all_h, all_z, all_P, ~] = run_isothermal_grading(comp, press_ref, temp, ...
    Pc, Tc, w, M*1e-3, BIP, h_ref, h_down, h_up, dh, R, g, tol, maxiter);

fprintf('  x_H2: %.4f (top) to %.4f (bot) mol%%\n', all_z(1,1)*100, all_z(end,1)*100);

figure('Position', [100 100 750 600]);
plot(all_z(:,1)*100, all_h, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8);
xlabel('H_2 (mol%)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse', 'FontSize', 11);
grid on;
title('Heavy NG — H_2 Compositional Grading (180 bar, 320.43 K)', 'FontSize', 12);


fprintf('\n=== HEAVY NG: GOC DETECTION ===\n');

[uGOC_depth, uGOC_pressure, uGOC_temp, uGOC_comp, uGOC_type, uGOC_info] = detect_GOC(...
    comp, press_ref, temp, h_ref, depth_range, ...
    Pc, Tc, w, BIP, M, ...
    0, Cp_dummy, Hig_dummy, vt_method, vt_params, ...
    'grading_model', 'isothermal', 'scan_step', 5, 'verbose', true);

fprintf('\n--- Summary ---\n');
if ~isnan(uGOC_depth)
    fprintf('  GOC: h = %.2f m, P = %.2f bar, type = %s\n', uGOC_depth, uGOC_pressure/1e5, uGOC_type);
else
    fprintf('  GOC: not found\n');
end

figure('Position', [100 100 800 600]);
hold on;

plot(uGOC_info.scan_P/1e5, uGOC_info.scan_depths, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8);

valid_sat = ~isnan(uGOC_info.scan_Psat);
if any(valid_sat)
    plot(uGOC_info.scan_Psat(valid_sat)/1e5, uGOC_info.scan_depths(valid_sat), '--', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end

if ~isnan(uGOC_depth)
    plot(uGOC_pressure/1e5, uGOC_depth, 'd', 'MarkerSize', 14, 'MarkerFaceColor', [0.49 0.18 0.56], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
end

xlabel('Pressure (bar)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse', 'FontSize', 11);
entries = {'P_{res}', 'P_{sat}'};
if ~isnan(uGOC_depth), entries{end+1} = 'Undersaturated GOC'; end
legend(entries, 'Location', 'best', 'FontSize', 10);
grid on;
title('Heavy NG — Reservoir Pressure vs Saturation Pressure (180 bar, 320.43 K)', 'FontSize', 12);

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================
function [all_h, all_z, all_P, all_rho] = run_isothermal_grading(comp_ref, press_ref, temp, Pc, Tc, w, M_kgmol, BIP, h_ref, h_down, h_up, dh, R, g, tol, maxiter)

    nc = length(comp_ref);
    n_down = length(h_down);
    n_up   = length(h_up);

    comp_down = zeros(n_down, nc); P_down = zeros(n_down,1); rho_down = zeros(n_down,1);
    comp_up   = zeros(n_up, nc);   P_up   = zeros(n_up,1);   rho_up   = zeros(n_up,1);

    [phi_ref, Z_ref] = fugacitycoef_multicomp_vapor(comp_ref, press_ref, temp, Pc, Tc, w, BIP, 'PR');
    rho_ref = press_ref * sum(comp_ref .* M_kgmol) / (Z_ref * R * temp);

    z_prev = comp_ref; P_prev = press_ref; phi_prev = phi_ref;
    for k = 1:n_down
        [~, Z_prev] = fugacitycoef_multicomp_vapor(z_prev, P_prev, temp, Pc, Tc, w, BIP, 'PR');
        rho_prev = P_prev * sum(z_prev .* M_kgmol) / (Z_prev * R * temp);
        P_guess  = P_prev + rho_prev * g * dh;
        fi_target = (phi_prev .* z_prev * P_prev) .* exp(M_kgmol * g * dh / (R * temp));
        z_new = z_prev; P_new = P_guess;
        for iter = 1:maxiter
            [phi_new, Z_new] = fugacitycoef_multicomp_vapor(z_new, P_new, temp, Pc, Tc, w, BIP, 'PR');
            z_calc = fi_target ./ (phi_new * P_new);
            z_calc = z_calc / sum(z_calc);
            rho_new = P_new * sum(z_calc .* M_kgmol) / (Z_new * R * temp);
            P_calc  = P_prev + rho_new * g * dh;
            if max(abs(z_calc - z_new)) < tol && abs(P_calc - P_new)/P_new < tol, break; end
            z_new = z_calc; P_new = P_calc;
        end
        comp_down(k,:) = z_new'; P_down(k) = P_new;
        rho_down(k) = P_new * sum(z_new .* M_kgmol) / (Z_new * R * temp);
        z_prev = z_new; P_prev = P_new; phi_prev = phi_new;
    end

    z_prev = comp_ref; P_prev = press_ref;
    [phi_prev, ~] = fugacitycoef_multicomp_vapor(comp_ref, press_ref, temp, Pc, Tc, w, BIP, 'PR');
    for k = 1:n_up
        [~, Z_prev] = fugacitycoef_multicomp_vapor(z_prev, P_prev, temp, Pc, Tc, w, BIP, 'PR');
        rho_prev = P_prev * sum(z_prev .* M_kgmol) / (Z_prev * R * temp);
        P_guess  = P_prev - rho_prev * g * dh;
        fi_target = (phi_prev .* z_prev * P_prev) .* exp(-M_kgmol * g * dh / (R * temp));
        z_new = z_prev; P_new = P_guess;
        for iter = 1:maxiter
            [phi_new, Z_new] = fugacitycoef_multicomp_vapor(z_new, P_new, temp, Pc, Tc, w, BIP, 'PR');
            z_calc = fi_target ./ (phi_new * P_new);
            z_calc = z_calc / sum(z_calc);
            rho_new = P_new * sum(z_calc .* M_kgmol) / (Z_new * R * temp);
            P_calc  = P_prev - rho_new * g * dh;
            if max(abs(z_calc - z_new)) < tol && abs(P_calc - P_new)/P_new < tol, break; end
            z_new = z_calc; P_new = P_calc;
        end
        comp_up(k,:) = z_new'; P_up(k) = P_new;
        rho_up(k) = P_new * sum(z_new .* M_kgmol) / (Z_new * R * temp);
        z_prev = z_new; P_prev = P_new; phi_prev = phi_new;
    end

    all_h   = [flipud(h_up); h_ref; h_down];
    all_z   = [flipud(comp_up); comp_ref'; comp_down];
    all_P   = [flipud(P_up); press_ref; P_down];
    all_rho = [flipud(rho_up); rho_ref; rho_down];
end