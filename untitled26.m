clear; clc;

R = 8.3144598;

P = 206 * 6894.757;
T = 288.75;

comp = [0.3094; 0.6906];
Tc = [305.322; 369.89];
Pc = [48.72e5; 42.48e5];
omega = [0.0995; 0.1521];
Zc = [0.2793; 0.2763];
MW = [30.069; 44.0956];
Vc = Zc .* R .* Tc ./ Pc * 1e6;
components = {'C2', 'C3'};
BIP = zeros(2);

V_exp = 83.573197;
rho_exp = 0.4757;
MW_mix = sum(comp .* MW);

vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

methods = [0 1 2 3 4 5 7 8 9 10];
names = {'PR (no VT)', 'PR + Peneloux', 'PR + Magoulas-Tassios', ...
         'PR + Ungerer-Batut', 'PR + Baled', 'PR + Abudour', ...
         'SRK (no VT)', 'SRK + Pina-Martinez', 'SRK + Chen-Li', ...
         'SRK + Baled'};

fprintf('==========================================================================\n');
fprintf('calculate_density VALIDATION: C2/C3 at 206 psi, T = 288.75 K\n');
fprintf('==========================================================================\n');
fprintf('Experimental: V = %.4f cm3/mol, rho = %.4f g/cm3 (%.2f kg/m3)\n\n', V_exp, rho_exp, rho_exp*1000);
fprintf('%-25s %6s %14s %12s %12s %10s\n', 'Method', 'VT#', 'V [cm3/mol]', 'rho [kg/m3]', 'rho [g/cm3]', 'Err [%]');
fprintf('%s\n', repmat('-', 1, 85));

for i = 1:length(methods)
    try
        [rho, V_m, Z, c_mix] = calculate_density(comp, P, T, Pc, Tc, omega, BIP, MW, methods(i), vt_params);
        V_cm3 = V_m * 1e6;
        rho_gcc = rho / 1000;
        err = (rho_gcc - rho_exp) / rho_exp * 100;
        fprintf('%-25s %6d %14.4f %12.2f %12.4f %10.2f\n', names{i}, methods(i), V_cm3, rho, rho_gcc, err);
    catch ME
        fprintf('%-25s %6d  *** ERROR: %s\n', names{i}, methods(i), ME.message);
    end
end

fprintf('%s\n', repmat('-', 1, 85));

