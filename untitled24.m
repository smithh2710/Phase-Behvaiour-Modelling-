clear; clc;

R = 8.314472;

P = 0.101352e6;
T = 20 + 273.15;

comp = 1;
Tc = 647.096;
Pc = 22.064e6;
omega = 0.3449;
Zc = 0.229;
MW = 18.015;
Vc = Zc * R * Tc / Pc * 1e6;
components = {'H2O'};
BIP = 0;

rho_exp = 998.29;
V_exp = MW / rho_exp * 1000;

vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

methods = [0 1 2 3 4 5 7 8 9 10];
names = {'PR (no VT)', 'PR + Peneloux', 'PR + Magoulas-Tassios', ...
         'PR + Ungerer-Batut', 'PR + Baled', 'PR + Abudour', ...
         'SRK (no VT)', 'SRK + Pina-Martinez', 'SRK + Chen-Li', ...
         'SRK + Baled'};

fprintf('==========================================================================\n');
fprintf('calculate_density VALIDATION: Water at 1 atm, T = 293.15 K\n');
fprintf('==========================================================================\n');
fprintf('Experimental: V = %.4f cm3/mol, rho = %.2f kg/m3\n\n', V_exp, rho_exp);
fprintf('%-25s %6s %14s %12s %10s\n', 'Method', 'VT#', 'V [cm3/mol]', 'rho [kg/m3]', 'Err [%]');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:length(methods)
    try
        [rho, V_m, Z, c_mix] = calculate_density(comp, P, T, Pc, Tc, omega, BIP, MW, methods(i), vt_params);
        V_cm3 = V_m * 1e6;
        err = (rho - rho_exp) / rho_exp * 100;
        fprintf('%-25s %6d %14.4f %12.2f %10.2f\n', names{i}, methods(i), V_cm3, rho, err);
    catch ME
        fprintf('%-25s %6d  *** ERROR: %s\n', names{i}, methods(i), ME.message);
    end
end


% [rho, V_m, Z, c_mix] = calculate_density(comp, P, T, Pc, Tc, omega, BIP, MW, 9, vt_params)
