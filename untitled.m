clear; clc; close all;

R = 8.3144598;

%% Component data: C2, C3
components = {'C2','C3'};
M_gmol = [30.069; 44.0956];
ncomp  = length(M_gmol);

%% SRK characterization
Tc_srk = [305.322; 369.89];
Pc_srk = [48.838e5; 42.46e5];
w_srk  = [0.0995; 0.1521];
BIP_srk = zeros(ncomp);

Zc_srk = 0.29056 - 0.08775 * w_srk;
vtp_srk.Zc = Zc_srk;
vtp_srk.Vc = Zc_srk .* R .* Tc_srk ./ Pc_srk * 1e6;
vtp_srk.components = components;

%% PR characterization
Tc_pr = [305.322; 369.89];
Pc_pr = [48.838e5; 42.46e5];
w_pr  = [0.0995; 0.1521];
BIP_pr = zeros(ncomp);

Zc_pr  = 0.29056 - 0.08775 * w_pr;
Vc_pr  = Zc_pr .* R .* Tc_pr ./ Pc_pr * 1e6;
ZRA_pr = Zc_pr;
c_peneloux_pr = 0.50033 * R * Tc_pr .* (0.25969 - ZRA_pr) ./ Pc_pr * 1e6;

vtp_pr.Zc = Zc_pr;
vtp_pr.Vc = Vc_pr;
vtp_pr.components = components;
vtp_pr.c_custom   = c_peneloux_pr;

%% Method definitions
all_methods = [0, 3, 4, 5, 6, 7, 8, 9];
all_names = {'PR', 'PR + Ungerer-Batut', ...
             'PR + Baled', 'PR + Abudour', 'PR + Constant VT', ...
             'SRK', 'SRK + Pina-Martinez', 'SRK + Chen-Li'};
nM = length(all_methods);
is_srk = all_methods >= 7;

%% Experimental data
% C2/C3 binary at T = 288.75 K (15.6 C), four pressures
T_exp     = 288.75;
P_exp_MPa = [1.42032000174; 1.651294370955; 1.8615844683; 2.919929712315];
x_C2      = [0.3094; 0.4037; 0.4855; 0.8635];
x_C3      = 1 - x_C2;
rho_exp   = [0.4757; 0.4649; 0.4530; 0.3897];   % g/cm3

nPts = length(P_exp_MPa);

%% Compute model densities
rho_model = NaN(nPts, nM);
for iPt = 1:nPts
    z = [x_C2(iPt); x_C3(iPt)];
    press = P_exp_MPa(iPt) * 1e6;
    Temp  = T_exp;
    for iM = 1:nM
        try
            if is_srk(iM)
                [rho_kg,~,~,~] = calculate_density(z, press, Temp, Pc_srk, Tc_srk, w_srk, BIP_srk, M_gmol, all_methods(iM), vtp_srk);
            else
                [rho_kg,~,~,~] = calculate_density(z, press, Temp, Pc_pr,  Tc_pr,  w_pr,  BIP_pr,  M_gmol, all_methods(iM), vtp_pr);
            end
            rho_model(iPt, iM) = rho_kg / 1000;
        catch
            rho_model(iPt, iM) = NaN;
        end
    end
end

%% Predicted density table
fprintf('\n=== Predicted density (g/cm3), C2/C3 binary, T = %.2f K ===\n', T_exp);
fprintf('%-25s', 'Method');
for iPt = 1:nPts
    fprintf('  P=%.2fMPa', P_exp_MPa(iPt));
end
fprintf('\n');
fprintf('%-25s', 'x_C2');
for iPt = 1:nPts
    fprintf('  x=%6.4f ', x_C2(iPt));
end
fprintf('\n');
fprintf('%-25s', 'rho_exp');
for iPt = 1:nPts
    fprintf('  %8.4f ', rho_exp(iPt));
end
fprintf('\n%s\n', repmat('-', 1, 25 + 11*nPts));
for iM = 1:nM
    fprintf('%-25s', all_names{iM});
    for iPt = 1:nPts
        if isnan(rho_model(iPt, iM))
            fprintf('  %8s ', '--');
        else
            fprintf('  %8.4f ', rho_model(iPt, iM));
        end
    end
    fprintf('\n');
end

%% Relative error per point + overall AARD + max error
relErr = abs(rho_model - rho_exp) ./ rho_exp * 100;

AARD_overall   = zeros(nM, 1);
maxErr_overall = zeros(nM, 1);
for iM = 1:nM
    valid = ~isnan(relErr(:, iM));
    if any(valid)
        AARD_overall(iM)   = mean(relErr(valid, iM));
        maxErr_overall(iM) = max(relErr(valid, iM));
    else
        AARD_overall(iM)   = NaN;
        maxErr_overall(iM) = NaN;
    end
end

fprintf('\n=== Relative error per point (%%), C2/C3 binary, T = %.2f K ===\n', T_exp);
fprintf('%-25s', 'Method');
for iPt = 1:nPts
    fprintf('  P=%.2fMPa', P_exp_MPa(iPt));
end
fprintf(' %10s %10s\n', 'AARD', 'Max');
fprintf('%s\n', repmat('-', 1, 25 + 11*nPts + 22));
for iM = 1:nM
    fprintf('%-25s', all_names{iM});
    for iPt = 1:nPts
        if isnan(relErr(iPt, iM))
            fprintf('  %8s ', '--');
        else
            fprintf('  %8.2f ', relErr(iPt, iM));
        end
    end
    fprintf(' %10.2f %10.2f\n', AARD_overall(iM), maxErr_overall(iM));
end

fprintf('\nN points: %d (all at T = %.2f K)\n', nPts, T_exp);