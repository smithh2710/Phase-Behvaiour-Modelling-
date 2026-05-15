function [H_mix_specific, H_partial_specific, H_ig_MR, H_res_MR, H_total_MR] = calculate_absolute_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components)

    R = 8.3144598;
    T_ref = 273.15;

    if nargin < 11 || isempty(eos_type)
        eos_type = 'PR';
    end
    if nargin < 12
        components = {};
    end

    n = length(comp);
    comp = comp(:);
    comp = comp / sum(comp);

    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    H_ig_ref = H_ig_ref(:);

    H_ig = zeros(n, 1);
    for i = 1:n
        C1 = Cp_coeffs(i, 1);
        C2 = Cp_coeffs(i, 2);
        C3 = Cp_coeffs(i, 3);
        C4 = Cp_coeffs(i, 4);
        
        delta_H_ig = C1 * (T - T_ref) + ...
                     C2 / 2 * (T^2 - T_ref^2) + ...
                     C3 / 3 * (T^3 - T_ref^3) + ...
                     C4 / 4 * (T^4 - T_ref^4);
        
        H_ig(i) = H_ig_ref(i) + delta_H_ig;
    end

    H_res = calculate_partial_molar_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R, eos_type, components);

    H_partial = H_ig + H_res;
    
    H_mix = sum(comp .* H_partial);
    
    M_mix = sum(comp .* M_gmol);
    
    H_mix_specific = H_mix / M_mix;
    
    H_partial_specific = H_partial ./ M_gmol;
    
    H_ig_MR = H_ig ./ (M_gmol * R);
    
    H_res_MR = H_res ./ (M_gmol * R);
    
    H_total_MR = H_partial ./ (M_gmol * R);

end


function H_res = calculate_partial_molar_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R, eos_type, components)

    n = length(comp);
    H_res = zeros(n, 1);

    dT = max(0.1, 1e-4 * T);

    [phi_minus, ~] = fugacitycoef_multicomp(comp, P, T - dT, Pc, Tc, acentric, BIP, eos_type, components);
    [phi_plus, ~]  = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP, eos_type, components);

    for i = 1:n
        if phi_plus(i) > 1e-15 && phi_minus(i) > 1e-15
            dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
            H_res(i) = -R * T^2 * dln_phi_dT;
        else
            [phi_0, ~] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type, components);
            [phi_fwd, ~] = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP, eos_type, components);
            if phi_0(i) > 1e-15 && phi_fwd(i) > 1e-15
                dln_phi_dT = (log(phi_fwd(i)) - log(phi_0(i))) / dT;
                H_res(i) = -R * T^2 * dln_phi_dT;
            else
                warning('Component %d has invalid fugacity coefficient', i);
                H_res(i) = 0;
            end
        end
    end

    bad_idx = ~isfinite(H_res);
    if any(bad_idx)
        warning('Non-finite residual enthalpies for %d components', sum(bad_idx));
        H_res(bad_idx) = 0;
    end
end