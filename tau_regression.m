clear; clc; close all;

%% ========================================================================
%  PEDERSEN 2024 - PER-COMPONENT TAU REGRESSION
%  Fix: N2=2.7, CO2=3.0, C1=2.5, C2=3.2, C3=3.5 (Ewell-Eyring)
%  Regress: iC4, nC4, iC5, nC5, C6, C7-C10, C11-C20, C21-C30, C31-C80
%  Bounds: [2.5, 6.0] for each
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

Cp_coeffs = [
    31.15,  -0.014,    2.68e-5,  -1.17e-8;
    19.79,   0.073,   -5.60e-5,   1.72e-8;
    19.25,   0.052,    1.20e-5,  -1.13e-8;
     5.41,   0.178,   -6.94e-5,   8.71e-9;
    -4.22,   0.306,   -1.59e-4,   3.21e-8;
    -1.39,   0.385,   -1.85e-4,   2.90e-8;
     9.49,   0.331,   -1.11e-4,  -2.82e-9;
     9.52,   0.507,   -2.73e-4,   5.72e-8;
    -3.63,   0.487,   -2.58e-4,   5.30e-8;
    -4.41,   0.582,   -3.12e-4,   6.49e-8;
     4.50,   0.620,   -2.20e-4,   0.000;
     7.50,   1.120,   -4.50e-4,   0.000;
    12.00,   1.920,   -7.70e-4,   0.000;
    20.00,   3.160,   -1.27e-3,   0.000];

H_ig_Jg = [0.0; 1231.4; 0.0; 304.3; 534.1; 620.3; 613.9; 461.5; 399.9; 393.3; ...
           147.1; 873.8; 883.2; 84.4];
H_ig_ref = H_ig_Jg .* M_gmol;

h_ref     = 2590;
press_ref = 284.5e5;
temp_ref  = 94.9 + 273.15;
dTdh      = 0.027;
tol = 1e-10;
maxiter = 500;

vt_method = 12;
vt_params = struct();
vt_params.c_custom = c_peneloux_cm3;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

idx_C1  = 3;
idx_C7p = 11:14;

%% ========================================================================
%  FIXED TAU VALUES (Ewell-Eyring literature)
%  ========================================================================
tau_fixed_idx = [1, 2, 3, 4, 5];                    % N2, CO2, C1, C2, C3
tau_fixed_val = [2.7; 3.0; 2.5; 3.2; 3.5];          % from Ewell-Eyring data

tau_free_idx  = [6, 7, 8, 9, 10, 11, 12, 13, 14];   % iC4 through C31-C80
n_free = length(tau_free_idx);

fprintf('=== PER-COMPONENT TAU REGRESSION ===\n');
fprintf('Fixed components:\n');
for i = 1:length(tau_fixed_idx)
    fprintf('  %-10s  tau = %.1f (Ewell-Eyring)\n', components{tau_fixed_idx(i)}, tau_fixed_val(i));
end
fprintf('Free components (%d): iC4, nC4, iC5, nC5, C6, C7-C10, C11-C20, C21-C30, C31-C80\n', n_free);
fprintf('Bounds: [2.5, 6.0]\n\n');

%% ========================================================================
%  EXPERIMENTAL DATA
%  ========================================================================
exp_h = [2431; 2559; 2590; 2618; 2635; 2636; 2658; 2660; 2697; 2712; 2728; ...
         2746; 2758; 2758; 2769];
exp_C1 = [75.747; 75.228; 49.946; 49.510; 49.387; 49.884; 48.890; 47.936; 47.234; 47.213; 47.096; ...
          45.028; 45.323; 45.322; 44.882];

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

exp_C7_C10 = [1.363+1.087+0.546+0.369; 1.374+1.176+0.557+0.378; ...
              4.006+3.969+2.303+1.822; 3.994+3.996+2.335+1.854; ...
              3.958+4.089+2.464+1.664; 3.719+3.840+2.358+1.854; ...
              3.851+3.970+2.408+1.940; 4.123+4.237+2.571+1.990; ...
              4.183+4.324+2.612+1.794; 4.096+4.207+2.605+1.996; ...
              4.069+4.215+2.593+2.016; 4.281+4.440+2.647+2.206; ...
              4.223+4.373+2.616+2.158; 4.190+4.260+2.589+2.184; ...
              4.229+4.387+2.630+2.197];

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

exp_light = [95.653; 95.463; 74.713; 74.276; 73.529; 74.360; 73.469; 72.636; 71.632; 72.069; 71.992; ...
             70.551; 70.933; 70.872; 70.348];
exp_C7p = 100 - exp_light;
exp_C11p = exp_C11_C19 + exp_C20p;

exp_Pres = [278.3; 282.4; 284.2; 285.8; 286.6; 286.7; 288.1; 288.2; 290.4; 291.3; 292.4; ...
            291.6; 292.4; 292.4; 293.1];
exp_Psat = [278.3; 282.4; 284.1; 280.0; 268.5; 267.4; 265.0; 259.2; 254.2; 255.3; 253.8; ...
            241.2; 241.6; 240.6; 239.1];
exp_GOR  = [3875; 3061; 278.9; 276.4; NaN; 280.9; 264.5; NaN; 235.5; 229.3; 228.7; ...
            220.7; 227.8; 221.7; 213.4];
is_gas = abs(exp_Pres - exp_Psat) < 2;

h_oil_targets = [2590, 2600, 2618, 2635, 2636, 2658, 2660, 2697, 2712, 2728, 2746, 2758, 2769];
h_gas_base = [2559, 2540, 2500, 2460, 2431];

%% ========================================================================
%  HELPER: Build full tau vector from free parameters
%  ========================================================================
build_tau = @(p_free) build_tau_vec(p_free, tau_fixed_idx, tau_fixed_val, tau_free_idx, n);

%% ========================================================================
%  STAGE 0: TAU=4 BASELINE
%  ========================================================================
fprintf('=== TAU=4 BASELINE ===\n');
tau_const = 4.0 * ones(n, 1);
res_tau4 = run_firoozabadi_forward(tau_const, comp_ref, h_ref, press_ref, temp_ref, ...
    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
    vt_method, vt_params, components, n, idx_C1, idx_C7p, ...
    h_oil_targets, h_gas_base, tol, maxiter);
fprintf('  GOC = %.1f m\n', res_tau4.GOC);

fobj_tau4 = eval_objective(res_tau4, exp_h, exp_C1, exp_C2C6, idx_C2C6, exp_C7_C10, exp_C11p, exp_C7p, is_gas);
fprintf('  Objective = %.6f\n\n', fobj_tau4);

%% ========================================================================
%  REGRESSION: 9 free tau values via fmincon
%  ========================================================================
fprintf('=== REGRESSION: 9 free tau values ===\n');

% Initial guess: tau=4 for all free components
p0 = 4.0 * ones(n_free, 1);
lb = 2.5 * ones(n_free, 1);
ub = 6.0 * ones(n_free, 1);

obj_handle = @(p_free) tau_obj_percomp(p_free, build_tau, ...
    comp_ref, h_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, ...
    M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, n, idx_C1, idx_C7p, ...
    h_oil_targets, exp_h, exp_C1, exp_C2C6, idx_C2C6, exp_C7_C10, exp_C11p, exp_C7p, is_gas);

fobj_init = obj_handle(p0);
fprintf('  Initial objective (all tau=4): %.6f\n', fobj_init);

opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 2000, ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-4, ...
    'FiniteDifferenceStepSize', 0.05, ...
    'Algorithm', 'interior-point');

[p_opt, fobj_opt, exitflag] = fmincon(obj_handle, p0, [], [], [], [], lb, ub, [], opts);

tau_opt = build_tau(p_opt);

fprintf('\n=== REGRESSION RESULTS ===\n');
fprintf('  Exit flag: %d\n', exitflag);
fprintf('  Objective: %.6f -> %.6f\n', fobj_init, fobj_opt);
fprintf('\n  %-10s  %8s  %8s  %8s  %8s\n', 'Component', 'M', 'tau=4', 'tau_opt', 'Status');
fprintf('  %s\n', repmat('-', 1, 55));
for i = 1:n
    if ismember(i, tau_fixed_idx)
        status = 'FIXED';
    else
        status = '';
        if abs(tau_opt(i) - 2.5) < 0.01, status = 'AT LB'; end
        if abs(tau_opt(i) - 6.0) < 0.01, status = 'AT UB'; end
    end
    fprintf('  %-10s  %8.1f  %8.3f  %8.3f  %s\n', components{i}, M_gmol(i), 4.0, tau_opt(i), status);
end

%% ========================================================================
%  FORWARD SIMULATION WITH OPTIMIZED TAU
%  ========================================================================
fprintf('\n=== Forward simulation with optimized tau ===\n');
res_opt = run_firoozabadi_forward(tau_opt, comp_ref, h_ref, press_ref, temp_ref, ...
    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
    vt_method, vt_params, components, n, idx_C1, idx_C7p, ...
    h_oil_targets, h_gas_base, tol, maxiter);
fprintf('  GOC = %.1f m\n', res_opt.GOC);

%% ========================================================================
%  PLOTS
%  ========================================================================
colors = {[0.49 0.18 0.56], [0.00 0.45 0.74]};
lstyles = {'-.', '-'};
lwidths = [1.5, 2.0];
names = {'Firoozabadi (\tau=4)', 'Firoozabadi (per-comp \tau)'};
all_res = {res_tau4, res_opt};

figure('Position', [50 50 1400 900]);

subplot(2,2,1); hold on;
for m = 1:2
    plot(all_res{m}.Psat/1e5, all_res{m}.h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(all_res{1}.P/1e5, all_res{1}.h, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(exp_Pres(is_gas), exp_h(is_gas), 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(exp_Psat(~is_gas), exp_h(~is_gas), 'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
plot(exp_Pres(~is_gas), exp_h(~is_gas), 'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('Pressure (bar)'); ylabel('Depth (m)');
legend([names, 'P_{res} gas', 'P_{sat} oil', 'P_{res} oil'], 'Location', 'best', 'FontSize', 7);
title('(a) P^{res} and P^{sat} vs. Depth');

subplot(2,2,2); hold on;
for m = 1:2
    valid = ~isnan(all_res{m}.GOR) & all_res{m}.GOR > 0;
    if any(valid)
        plot(all_res{m}.GOR(valid), all_res{m}.h(valid), 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
    end
end
valid_gor = ~isnan(exp_GOR);
plot(exp_GOR(is_gas & valid_gor), exp_h(is_gas & valid_gor), 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(exp_GOR(~is_gas & valid_gor), exp_h(~is_gas & valid_gor), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('GOR (Sm^3/Sm^3)'); ylabel('Depth (m)');
legend([names, 'Gas exp', 'Oil exp'], 'Location', 'best', 'FontSize', 7);
title('(b) GOR vs. Depth');

subplot(2,2,3); hold on;
for m = 1:2
    plot(all_res{m}.C1, all_res{m}.h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C1(is_gas), exp_h(is_gas), 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(exp_C1(~is_gas), exp_h(~is_gas), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('C_1 (mol%)'); ylabel('Depth (m)');
legend([names, 'Gas', 'Oil'], 'Location', 'best');
title('(c) C_1 mol% vs. Depth');

subplot(2,2,4); hold on;
for m = 1:2
    plot(all_res{m}.C7p, all_res{m}.h, 'Color', colors{m}, 'LineStyle', lstyles{m}, 'LineWidth', lwidths(m));
end
plot(exp_C7p(is_gas), exp_h(is_gas), 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(exp_C7p(~is_gas), exp_h(~is_gas), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
set(gca, 'YDir', 'reverse'); grid on;
xlabel('C_{7+} (mol%)'); ylabel('Depth (m)');
legend([names, 'Gas zone', 'Oil zone'], 'Location', 'best');
title('(d) C_{7+} mol% vs. Depth');

sgtitle('Per-component \tau regression (5 fixed, 9 free)', 'FontSize', 13);

%% ========================================================================
%  RMSE TABLE
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  GOC: tau=4 = %.1f m  |  per-comp = %.1f m  |  exp: 2559-2590 m\n', res_tau4.GOC, res_opt.GOC);
fprintf('==========================================================\n');
fprintf('  %-35s  %8s  %8s  %9s  %8s  %9s\n', 'Model', 'RMSE(C1)', 'RMSE(C7+)', 'RMSE(C7-10)', 'RMSE(C11+)', 'RMSE(C2-C6)');
fprintf('  %s\n', repmat('-', 1, 95));
all_labels = {'Firoozabadi (tau=4)', 'Firoozabadi (per-comp tau)'};
all_r = {res_tau4, res_opt};
for m = 1:2
    err_C1 = []; err_C7p = []; err_C710 = []; err_C11p = []; err_C2C6_all = [];
    for k = 1:length(exp_h)
        if is_gas(k), continue; end
        [~, idx] = min(abs(all_r{m}.h - exp_h(k)));
        if abs(all_r{m}.h(idx) - exp_h(k)) < 15
            err_C1   = [err_C1;   all_r{m}.C1(idx) - exp_C1(k)];
            err_C7p  = [err_C7p;  all_r{m}.C7p(idx) - exp_C7p(k)];
            err_C710 = [err_C710; all_r{m}.comp(idx,11)*100 - exp_C7_C10(k)];
            err_C11p = [err_C11p; sum(all_r{m}.comp(idx,12:14))*100 - exp_C11p(k)];
            for j = 1:length(idx_C2C6)
                err_C2C6_all = [err_C2C6_all; all_r{m}.comp(idx,idx_C2C6(j))*100 - exp_C2C6(k,j)];
            end
        end
    end
    if ~isempty(err_C1)
        fprintf('  %-35s  %7.2f %%  %7.2f %%  %9.2f %%  %8.2f %%  %9.2f %%\n', all_labels{m}, ...
            sqrt(mean(err_C1.^2)), sqrt(mean(err_C7p.^2)), sqrt(mean(err_C710.^2)), ...
            sqrt(mean(err_C11p.^2)), sqrt(mean(err_C2C6_all.^2)));
    end
end

%% ========================================================================
%  TAU BAR CHART
%  ========================================================================
figure('Position', [200 200 800 400]);
bar_data = [4*ones(n,1), tau_opt];
b = bar(bar_data);
b(1).FaceColor = [0.49 0.18 0.56];
b(2).FaceColor = [0.00 0.45 0.74];
set(gca, 'XTickLabel', components, 'XTickLabelRotation', 45);
ylabel('\tau'); grid on;
legend('\tau = 4 (constant)', '\tau_{opt} (per-component)', 'Location', 'northwest');
title('Per-component \tau values: constant vs optimized');
yline(2.5, 'r--', 'Lower bound', 'LabelHorizontalAlignment', 'left');
yline(6.0, 'r--', 'Upper bound', 'LabelHorizontalAlignment', 'left');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function tau_vec = build_tau_vec(p_free, fixed_idx, fixed_val, free_idx, n)
    tau_vec = zeros(n, 1);
    for i = 1:length(fixed_idx)
        tau_vec(fixed_idx(i)) = fixed_val(i);
    end
    for i = 1:length(free_idx)
        tau_vec(free_idx(i)) = p_free(i);
    end
end

function fobj = eval_objective(res, exp_h, exp_C1, exp_C2C6, idx_C2C6, exp_C7_C10, exp_C11p, exp_C7p, is_gas)
    fobj = 0;
    idx_C710  = 11;
    idx_C11p  = 12:14;
    for k = 1:length(exp_h)
        if is_gas(k), continue; end
        [~, idx] = min(abs(res.h - exp_h(k)));
        if abs(res.h(idx) - exp_h(k)) > 15, continue; end

        if exp_C1(k) > 0.1
            fobj = fobj + ((res.C1(idx) - exp_C1(k)) / exp_C1(k))^2;
        end

        for j = 1:length(idx_C2C6)
            sim_j = res.comp(idx, idx_C2C6(j)) * 100;
            exp_j = exp_C2C6(k, j);
            if exp_j > 0.05
                fobj = fobj + ((sim_j - exp_j) / exp_j)^2;
            end
        end

        sim_C710 = res.comp(idx, idx_C710) * 100;
        if exp_C7_C10(k) > 0.1
            fobj = fobj + 3 * ((sim_C710 - exp_C7_C10(k)) / exp_C7_C10(k))^2;
        end

        sim_C11p = sum(res.comp(idx, idx_C11p)) * 100;
        if exp_C11p(k) > 0.1
            fobj = fobj + 3 * ((sim_C11p - exp_C11p(k)) / exp_C11p(k))^2;
        end
    end
end

function fobj = tau_obj_percomp(p_free, build_tau_fn, ...
    comp_ref, h_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, ...
    M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, n, idx_C1, idx_C7p, ...
    h_oil_targets, exp_h, exp_C1, exp_C2C6, idx_C2C6, exp_C7_C10, exp_C11p, exp_C7p, is_gas)

    tau_vec = build_tau_fn(p_free);
    try
        res = run_firoozabadi_oil_only(tau_vec, comp_ref, h_ref, press_ref, temp_ref, ...
            dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
            vt_method, vt_params, n, idx_C1, idx_C7p, h_oil_targets);
        fobj = eval_objective(res, exp_h, exp_C1, exp_C2C6, idx_C2C6, exp_C7_C10, exp_C11p, exp_C7p, is_gas);
    catch
        fobj = 1e6;
    end
end

function res = run_firoozabadi_oil_only(tau_vec, comp_ref, h_ref, press_ref, temp_ref, ...
    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
    vt_method, vt_params, n, idx_C1, idx_C7p, h_oil_targets)

    n_oil = length(h_oil_targets);
    oil_comp = zeros(n_oil, n);

    for i = 1:n_oil
        [ch, ~, ~, ~, ~] = main_firoozabadi(h_oil_targets(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_vec, vt_method, vt_params);
        oil_comp(i,:) = ch';
    end

    res.h    = h_oil_targets(:);
    res.C1   = oil_comp(:, idx_C1) * 100;
    res.C7p  = sum(oil_comp(:, idx_C7p), 2) * 100;
    res.comp = oil_comp;
end

function res = run_firoozabadi_forward(tau_vec, comp_ref, h_ref, press_ref, temp_ref, ...
    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
    vt_method, vt_params, components, n, idx_C1, idx_C7p, ...
    h_oil_targets, h_gas_base, tol, maxiter)

    [GOC_h, ~, ~, GOC_T] = detect_sgoc(3, comp_ref, press_ref, temp_ref, h_ref, dTdh, ...
        [2500 2650], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau_vec(1));

    h_oil = [GOC_h, GOC_h+2, h_oil_targets];
    n_oil = length(h_oil);
    oil_P = zeros(n_oil,1); oil_Pbub = zeros(n_oil,1);
    oil_comp = zeros(n_oil, n); oil_rho = zeros(n_oil,1);
    oil_GOR = NaN(n_oil,1);

    for i = 1:n_oil
        [ch,Ph,Th,~,~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_vec, vt_method, vt_params);
        oil_P(i) = Ph; oil_comp(i,:) = ch';
        try
            [Pb_ext, ~] = pressbub_multicomp_newton(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
            oil_Pbub(i) = Pb_ext;
        catch, oil_Pbub(i) = NaN; end
        try, oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
        catch, oil_rho(i) = NaN; end
        try, oil_GOR(i) = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, ...
                'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'SRK');
        catch, oil_GOR(i) = NaN; end
    end

    [~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, ...
        Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);

    h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
    n_gas = length(h_gas);
    gas_P = zeros(n_gas,1); gas_Pdew = zeros(n_gas,1);
    gas_comp = zeros(n_gas, n); gas_rho = zeros(n_gas,1);
    gas_GOR = NaN(n_gas,1);

    for i = 1:n_gas
        [ch,Ph,Th,~,~] = main_firoozabadi(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau_vec, vt_method, vt_params);
        gas_P(i) = Ph; gas_comp(i,:) = ch';
        try
            [Pd_ext, ~] = pressdew_multicomp_ss(ch, Ph, Th, Pc, Tc, acentric, BIP, tol, maxiter, 'SRK', components);
            gas_Pdew(i) = Pd_ext;
        catch, gas_Pdew(i) = NaN; end
        try, gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
        catch, gas_rho(i) = NaN; end
        try, gas_GOR(i) = calculate_GOR_STO(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, ...
                'vt_method', vt_method, 'vt_params', vt_params, 'eos_type', 'SRK');
        catch, gas_GOR(i) = NaN; end
    end

    res.GOC  = GOC_h;
    res.h    = [flipud(h_gas'); h_oil'];
    res.P    = [flipud(gas_P); oil_P];
    res.Psat = [flipud(gas_Pdew); oil_Pbub];
    res.rho  = [flipud(gas_rho); oil_rho];
    res.C1   = [flipud(gas_comp(:,idx_C1)); oil_comp(:,idx_C1)] * 100;
    res.C7p  = [flipud(sum(gas_comp(:,idx_C7p),2)); sum(oil_comp(:,idx_C7p),2)] * 100;
    res.comp = [flipud(gas_comp); oil_comp];
    res.GOR  = [flipud(gas_GOR); oil_GOR];
end