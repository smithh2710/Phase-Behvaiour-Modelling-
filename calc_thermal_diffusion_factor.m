function [alpha_T, Q_star, dmu_dx] = calc_thermal_diffusion_factor(T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components, model, varargin)
%
% Thermal diffusion factor (Soret coefficient) for multicomponent mixtures.
%
% Definition (Firoozabadi, 2016):
%
%   alpha_T(i,j) = (Q_j^{*,M} - Q_i^{*,M}) / (x_i * dmu_i/dx_i)
%
% where Q_i^{*,M} is the net molar heat of transport from the selected
% model, and dmu_i/dx_i is the chemical potential derivative at constant
% T, P, and all other x_k.
%
% INPUTS:
%   T, P          - Temperature [K], Pressure [Pa]
%   comp          - Mole fractions [n x 1]
%   Pc, Tc        - Critical properties [Pa], [K]
%   acentric      - Acentric factors
%   BIP           - Binary interaction parameter matrix
%   M_gmol        - Molecular weights [g/mol]
%   Cp_coeffs     - Ideal gas Cp polynomial coefficients [n x 4]
%   H_ig_ref      - Ideal gas reference enthalpies at 273.15 K [J/mol]
%   eos_type      - 'PR' or 'SRK'
%   components    - Cell array of component names
%   model         - 'haase', 'kempers', or 'firoozabadi'
%
% Optional name-value pairs:
%   'tau'         - Firoozabadi tau parameter (default: 4.0)
%   'vt_method'   - Volume translation method (default: 0 = none)
%   'vt_params'   - VT parameter struct
%   'pair'        - [i, j] indices for specific binary pair (default: all)
%   'dx'          - Perturbation for numerical derivatives (default: 1e-5)
%
% OUTPUTS:
%   alpha_T  - Thermal diffusion factor matrix [n x n]
%              alpha_T(i,j) = (Q_j - Q_i) / (x_i * dmu_i/dx_i)
%   Q_star   - Net heat of transport for each component [J/mol], [n x 1]
%   dmu_dx   - Chemical potential derivatives [J/mol], [n x 1]
%              dmu_dx(i) = x_i * (d mu_i / d x_i)_{T,P,x_{k!=i}}

p = inputParser;
addParameter(p, 'tau', 4.0);
addParameter(p, 'vt_method', 0);
addParameter(p, 'vt_params', struct());
addParameter(p, 'pair', []);
addParameter(p, 'dx', 1e-5);
parse(p, varargin{:});
tau = p.Results.tau;
vt_method = p.Results.vt_method;
vt_params = p.Results.vt_params;
dx = p.Results.dx;

R = 8.3144598;
n = length(comp);
comp = comp(:); comp = comp / sum(comp);

%% ===== Q_i^{*,M}: Net heat of transport =====
switch lower(model)

    case 'haase'
        % Q_i^{*,Haase} = (M_i/M)*H_mix - H_tilde_i
        % where H_tilde_i = partial molar enthalpy, H_mix = molar enthalpy
        [H_mix_spec, H_partial_spec] = calculate_absolute_enthalpy(...
            T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components);

        M_mix = sum(comp .* M_gmol);
        H_mix = H_mix_spec * M_mix;           % J/mol
        H_tilde = H_partial_spec .* M_gmol;   % J/mol (partial molar)

        Q_star = (M_gmol / M_mix) .* H_mix - H_tilde;

    case 'kempers'
        % Q_i^{*,Kempers} = (V_bar_i / V_m)*H_mix - H_tilde_i
        [H_mix_spec, H_partial_spec] = calculate_absolute_enthalpy(...
            T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components);

        M_mix = sum(comp .* M_gmol);
        H_mix = H_mix_spec * M_mix;
        H_tilde = H_partial_spec .* M_gmol;

        [~, V_m] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type, components);
        V_m = V_m * R * T / P;  % if fugacitycoef returns Z

        V_partial = calc_partial_molar_volume(T, P, comp, Pc, Tc, acentric, BIP, eos_type, components, dx);

        Q_star = (V_partial / V_m) .* H_mix - H_tilde;

    case 'firoozabadi'
        % Q_i^{*,Firooz} = (1/tau) * (-delta_U_i + delta_U_mix * V_bar_i / V_m)
        % delta_U_i = H_res_i - P*V_bar_i  (residual internal energy)
        [~, ~, ~, H_res_MR] = calculate_absolute_enthalpy(...
            T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components);

        H_res = H_res_MR .* M_gmol * R;  % J/mol

        [~, Z] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type, components);
        V_m = Z * R * T / P;

        V_partial = calc_partial_molar_volume(T, P, comp, Pc, Tc, acentric, BIP, eos_type, components, dx);

        delta_U = H_res - P * V_partial;
        delta_U_mix = sum(comp .* delta_U);

        Q_star = (1/tau) * (-delta_U + delta_U_mix * V_partial / V_m);

    otherwise
        error('Unknown model: %s. Use haase, kempers, or firoozabadi.', model);
end

%% ===== dmu_i/dx_i from EOS =====
% mu_i = mu_i^0(T) + RT ln(phi_i * x_i * P)
% d(mu_i)/d(x_i) = RT * (1/x_i + d ln(phi_i)/d x_i)
% We compute x_i * d(mu_i)/d(x_i) = RT * (1 + x_i * d ln(phi_i)/d x_i)

[phi_0, ~] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type, components);
ln_phi_0 = log(phi_0(:));

dmu_dx = zeros(n, 1);

for i = 1:n
    comp_pert = comp;
    h = max(dx, dx * comp(i));

    if comp(i) + h > 1 - 1e-10
        h = -h;
    end

    comp_pert(i) = comp(i) + h;

    % Redistribute perturbation to maintain sum = 1
    others = true(n,1); others(i) = false;
    sum_others = sum(comp(others));
    if sum_others > 1e-15
        comp_pert(others) = comp(others) * (1 - comp_pert(i)) / sum_others;
    end
    comp_pert = max(comp_pert, 1e-30);
    comp_pert = comp_pert / sum(comp_pert);

    [phi_pert, ~] = fugacitycoef_multicomp(comp_pert, P, T, Pc, Tc, acentric, BIP, eos_type, components);
    ln_phi_pert = log(phi_pert(:));

    dlnphi_dxi = (ln_phi_pert(i) - ln_phi_0(i)) / h;

    % x_i * (dmu_i/dx_i) = RT * (1 + x_i * dlnphi/dxi)
    dmu_dx(i) = R * T * (1 + comp(i) * dlnphi_dxi);
end

%% ===== Thermal diffusion factor matrix =====
alpha_T_max = 100;
alpha_T = zeros(n, n);
for i = 1:n
    for j = 1:n
        if abs(dmu_dx(i)) > 1e-10
            raw = (Q_star(j) - Q_star(i)) / dmu_dx(i);
            alpha_T(i,j) = max(-alpha_T_max, min(alpha_T_max, raw));
        else
            alpha_T(i,j) = sign(Q_star(j) - Q_star(i)) * alpha_T_max;
        end
    end
end

end


%% ===== PARTIAL MOLAR VOLUME (numerical) =====
function V_partial = calc_partial_molar_volume(T, P, comp, Pc, Tc, acentric, BIP, eos_type, components, dx)

R = 8.3144598;
n = length(comp);
V_partial = zeros(n, 1);

[~, Z0] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type, components);
V_m = Z0 * R * T / P;

for i = 1:n
    % Add dn_i moles to system, compute dV/dn_i at constant T, P, n_{j!=i}
    dn = max(dx, dx * comp(i));
    n_total = 1.0;

    n_moles = comp * n_total;
    n_moles_pert = n_moles;
    n_moles_pert(i) = n_moles(i) + dn;
    n_total_pert = sum(n_moles_pert);
    comp_pert = n_moles_pert / n_total_pert;

    [~, Z_pert] = fugacitycoef_multicomp(comp_pert, P, T, Pc, Tc, acentric, BIP, eos_type, components);
    V_total_pert = n_total_pert * Z_pert * R * T / P;
    V_total_0 = n_total * V_m;

    V_partial(i) = (V_total_pert - V_total_0) / dn;
end

end