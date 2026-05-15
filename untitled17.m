%% VT Sign Convention Validation
%  Pure-component density check
%  Methane:  110 K, 50 bar    rho_exp = 0.4347 g/cm3 (NIST, compressed liquid)
%  n-Octane: 298.15 K, 1 atm  rho_exp = 0.6986 g/cm3 (DIPPR)
%
%  Both conditions chosen so the PR cubic has only one real root
%  (avoids vapor/liquid root selection ambiguity).
%
%  Expected behavior if sign convention is correct:
%    Methane (light, Tr=0.58):
%      no-VT PR over-predicts density --> Abudour should DECREASE density
%    n-Octane (heavy, Tr=0.52):
%      no-VT PR under-predicts density --> Abudour should INCREASE density
%
%  If both move toward experimental, sign convention is correct.

clear; clc;

R = 8.3144598;

%% ===== Test 1: Methane at 110 K, 50 bar (compressed liquid) =====
fprintf('\n========== TEST 1: PURE METHANE ==========\n');
fprintf('Conditions: T = 110 K, P = 50 bar (compressed liquid)\n');

T1   = 110;
P1   = 50e5;
comp = 1;
M    = 16.043;

Pc1 = 46.00e5;
Tc1 = 190.56;
w1  = 0.0115;
Zc1 = 0.29056 - 0.08775 * w1;
Vc1 = Zc1 * R * Tc1 / Pc1 * 1e6;

components1 = {'C1'};
BIP1        = 0;

vt_params           = struct();
vt_params.Vc        = Vc1;
vt_params.Zc        = Zc1;
vt_params.components = components1;

rho_exp_C1 = 0.4347;
V_exp_C1   = M / rho_exp_C1;

[rho_noVT, ~, ~, ~] = calculate_density(comp, P1, T1, Pc1, Tc1, w1, BIP1, M, 0, vt_params);
rho_noVT_C1 = rho_noVT / 1000;
V_noVT_C1   = M / rho_noVT_C1;

[rho_abu, ~, ~, ~]  = calculate_density(comp, P1, T1, Pc1, Tc1, w1, BIP1, M, 5, vt_params);
rho_abu_C1  = rho_abu / 1000;
V_abu_C1    = M / rho_abu_C1;

err_noVT_C1 = (rho_noVT_C1 - rho_exp_C1) / rho_exp_C1 * 100;
err_abu_C1  = (rho_abu_C1  - rho_exp_C1) / rho_exp_C1 * 100;

fprintf('\n%-20s %12s %12s %12s\n', 'Method', 'V (cm3/mol)', 'rho (g/cm3)', 'Error (%)');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%-20s %12.3f %12.4f %12s\n', 'Experimental', V_exp_C1, rho_exp_C1, '---');
fprintf('%-20s %12.3f %12.4f %+12.2f\n', 'PR (no VT)',  V_noVT_C1, rho_noVT_C1, err_noVT_C1);
fprintf('%-20s %12.3f %12.4f %+12.2f\n', 'PR + Abudour', V_abu_C1, rho_abu_C1, err_abu_C1);

if abs(err_abu_C1) < abs(err_noVT_C1)
    fprintf('\n  --> Abudour MOVED TOWARD experiment for C1.\n');
else
    fprintf('\n  --> Abudour MOVED AWAY FROM experiment for C1. Possible sign issue.\n');
end

%% ===== Test 2: n-Octane at 25 C, 1 atm =====
fprintf('\n\n========== TEST 2: PURE n-OCTANE ==========\n');
fprintf('Conditions: T = 298.15 K, P = 1.01325 bar (subcooled liquid)\n');

T2   = 298.15;
P2   = 1.01325e5;
comp = 1;
M    = 114.231;

Pc2 = 24.86e5;
Tc2 = 568.70;
w2  = 0.3996;
Zc2 = 0.29056 - 0.08775 * w2;
Vc2 = Zc2 * R * Tc2 / Pc2 * 1e6;

components2 = {'nC8'};
BIP2        = 0;

vt_params           = struct();
vt_params.Vc        = Vc2;
vt_params.Zc        = Zc2;
vt_params.components = components2;

rho_exp_C8 = 0.6986;
V_exp_C8   = M / rho_exp_C8;

[rho_noVT, ~, ~, ~] = calculate_density(comp, P2, T2, Pc2, Tc2, w2, BIP2, M, 0, vt_params);
rho_noVT_C8 = rho_noVT / 1000;
V_noVT_C8   = M / rho_noVT_C8;

[rho_abu, ~, ~, ~]  = calculate_density(comp, P2, T2, Pc2, Tc2, w2, BIP2, M, 5, vt_params);
rho_abu_C8  = rho_abu / 1000;
V_abu_C8    = M / rho_abu_C8;

err_noVT_C8 = (rho_noVT_C8 - rho_exp_C8) / rho_exp_C8 * 100;
err_abu_C8  = (rho_abu_C8  - rho_exp_C8) / rho_exp_C8 * 100;

fprintf('\n%-20s %12s %12s %12s\n', 'Method', 'V (cm3/mol)', 'rho (g/cm3)', 'Error (%)');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%-20s %12.3f %12.4f %12s\n', 'Experimental', V_exp_C8, rho_exp_C8, '---');
fprintf('%-20s %12.3f %12.4f %+12.2f\n', 'PR (no VT)',  V_noVT_C8, rho_noVT_C8, err_noVT_C8);
fprintf('%-20s %12.3f %12.4f %+12.2f\n', 'PR + Abudour', V_abu_C8, rho_abu_C8, err_abu_C8);

if abs(err_abu_C8) < abs(err_noVT_C8)
    fprintf('\n  --> Abudour MOVED TOWARD experiment for nC8.\n');
else
    fprintf('\n  --> Abudour MOVED AWAY FROM experiment for nC8. Possible sign issue.\n');
end

%% ===== Diagnosis =====
fprintf('\n\n========== DIAGNOSIS ==========\n');
C1_better  = abs(err_abu_C1) < abs(err_noVT_C1);
C8_better  = abs(err_abu_C8) < abs(err_noVT_C8);

if C1_better && C8_better
    fprintf('Both improved. Sign convention is correct.\n');
    fprintf('Abudour reducing liquid density for the Schulte fluid is therefore\n');
    fprintf('a real physical result, not a code bug.\n');
elseif ~C1_better && ~C8_better
    fprintf('Both got worse. Likely sign error in V_corrected = V_EOS - c.\n');
    fprintf('Check the sign of c in calculate_density.m.\n');
elseif ~C1_better && C8_better
    fprintf('Heavy works, light does not. Possible parameter miscalibration in\n');
    fprintf('Abudour, or distance function instability near critical conditions.\n');
elseif C1_better && ~C8_better
    fprintf('Light works, heavy does not. Inverse of typical bug. Check parameters.\n');
end

fprintf('\nSign of no-VT errors tells you PR''s native bias direction:\n');
if err_noVT_C1 > 0
    fprintf('  Methane: PR over-predicts density (positive error).\n');
else
    fprintf('  Methane: PR under-predicts density (negative error).\n');
end
if err_noVT_C8 > 0
    fprintf('  n-Octane: PR over-predicts density (positive error).\n');
else
    fprintf('  n-Octane: PR under-predicts density (negative error).\n');
end
fprintf('A correct VT should reduce |error| in both cases.\n');