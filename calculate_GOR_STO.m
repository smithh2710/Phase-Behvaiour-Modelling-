function [GOR, Bo, Bg, info] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, acentric, BIP, M_gmol, varargin)

    %% Parse inputs
    p = inputParser;
    p.FunctionName = 'calculate_GOR_STO';
    
    addRequired(p, 'comp', @isnumeric);
    addRequired(p, 'P_res', @isnumeric);
    addRequired(p, 'T_res', @isnumeric);
    addRequired(p, 'Pc', @isnumeric);
    addRequired(p, 'Tc', @isnumeric);
    addRequired(p, 'acentric', @isnumeric);
    addRequired(p, 'BIP', @isnumeric);
    addRequired(p, 'M_gmol', @isnumeric);
    
    addParameter(p, 'vt_method', 0);
    addParameter(p, 'vt_params', struct());
    addParameter(p, 'eos_type', 'PR', @ischar);
    addParameter(p, 'T_STO', 288.15, @isnumeric);
    addParameter(p, 'P_STO', 101325, @isnumeric);
    addParameter(p, 'verbose', false, @islogical);
    
    parse(p, comp, P_res, T_res, Pc, Tc, acentric, BIP, M_gmol, varargin{:});
    opts = p.Results;
    
    vt_method = opts.vt_method;
    vt_params = opts.vt_params;
    eos_type = upper(opts.eos_type);
    T_STO = opts.T_STO;
    P_STO = opts.P_STO;
    verbose = opts.verbose;
    
    %% Validate EOS type
    if ~ismember(eos_type, {'PR', 'SRK'})
        error('eos_type must be ''PR'' or ''SRK''. Got: %s', eos_type);
    end
    
    %% Initialize
    R = 8.3144598;
    
    comp = comp(:);
    n = length(comp);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    
    comp = comp / sum(comp);
    
    GOR = NaN;
    Bo = NaN;
    Rs = NaN;
    Bg = NaN;
    
    info = struct();
    info.eos_type = eos_type;
    info.P_STO = P_STO;
    info.T_STO = T_STO;
    info.P_res = P_res;
    info.T_res = T_res;
    
    
    %% Step 1: Flash to Stock Tank Conditions
    try
        [K_sto, comp_vap_sto, comp_liq_sto, V_frac] = vaporliquideq(P_STO, T_STO,comp, Pc, Tc, acentric, BIP, 1e-10, 200, eos_type);
        
        flash_converged = true;
        
        if V_frac < 1e-10
            V_frac = 0;
            comp_liq_sto = comp;
            comp_vap_sto = comp;
        elseif V_frac > (1 - 1e-10)
            if verbose
                fprintf('EOS flash returned all vapor (V=%.6f). Retrying with Wilson K-values...\n', V_frac);
            end
            K_wilson = wilsoneq(P_STO, T_STO, Pc, Tc, acentric);
            V_frac = solve_RR(comp, K_wilson);
            if V_frac > (1 - 1e-10)
                V_frac = 1 - 1e-10;
            end
            comp_vap_sto = K_wilson .* comp ./ (1 + V_frac * (K_wilson - 1));
            comp_liq_sto = comp ./ (1 + V_frac * (K_wilson - 1));
            if verbose
                fprintf('Wilson fallback: V_frac = %.6f, n_liq = %.6e\n', V_frac, 1-V_frac);
            end
        end
        
    catch ME
        flash_converged = false;
        if verbose
            fprintf('Flash failed: %s\n', ME.message);
            fprintf('Using Wilson K-values as fallback...\n');
        end
        
        K_wilson = wilsoneq(P_STO, T_STO, Pc, Tc, acentric);
        V_frac = solve_RR(comp, K_wilson);
        
        if V_frac < 0
            V_frac = 0;
        elseif V_frac > (1 - 1e-10)
            V_frac = 1 - 1e-10;
        end
        
        comp_vap_sto = K_wilson .* comp ./ (1 + V_frac * (K_wilson - 1));
        comp_liq_sto = comp ./ (1 + V_frac * (K_wilson - 1));
    end
    
    comp_vap_sto = comp_vap_sto(:);
    comp_liq_sto = comp_liq_sto(:);
    comp_vap_sto = comp_vap_sto / sum(comp_vap_sto);
    comp_liq_sto = comp_liq_sto / sum(comp_liq_sto);
    
    info.V_frac = V_frac;
    info.comp_vap_sto = comp_vap_sto;
    info.comp_liq_sto = comp_liq_sto;
    info.flash_converged = flash_converged;
    
    if verbose
        fprintf('Flash Results:\n');
        fprintf('  Vapor fraction: %.4f\n', V_frac);
        fprintf('  Liquid fraction: %.4f\n', 1 - V_frac);
        if V_frac > 0 && V_frac < 1
            fprintf('  Vapor C1: %.2f mol%%\n', comp_vap_sto(3)*100);
            fprintf('  Liquid C1: %.2f mol%%\n', comp_liq_sto(3)*100);
        end
        fprintf('\n');
    end
    
    %% Step 2: Calculate Molar Volumes at STO
    
    V_gas_molar = R * T_STO / P_STO;
    
    try
        [~, Z_liq] = fugacitycoef_multicomp_liquid(comp_liq_sto, P_STO, T_STO,Pc, Tc, acentric, BIP, eos_type);
        V_liq_molar_EOS = Z_liq * R * T_STO / P_STO;
        
        c_VT = calc_volume_translation(comp_liq_sto, P_STO, T_STO, Pc, Tc, acentric, vt_method, vt_params, n, eos_type, M_gmol,BIP);
        c_mix = sum(comp_liq_sto .* c_VT);
        V_liq_molar = V_liq_molar_EOS - c_mix;
        
        if V_liq_molar <= 0
            V_liq_molar = V_liq_molar_EOS;
        end
        
    catch ME
        if verbose
            fprintf('EOS liquid volume failed: %s\n', ME.message);
            fprintf('Using density estimate...\n');
        end
        
        MW_liq = sum(comp_liq_sto .* M_gmol);
        rho_liq = 700;
        V_liq_molar = MW_liq / (rho_liq * 1000);
    end
    
    info.V_gas_molar = V_gas_molar;
    info.V_liq_molar = V_liq_molar;
    
    
    %% Step 3: Calculate GOR
    
    n_gas = V_frac;
    n_liq = 1 - V_frac;
    
    V_gas_std = n_gas * V_gas_molar;
    V_liq_std = n_liq * V_liq_molar;
    
    if V_liq_std > 1e-15
        GOR = V_gas_std / V_liq_std;
    else
        GOR = Inf;
    end
    
    info.n_gas = n_gas;
    info.n_liq = n_liq;
    info.V_gas_std = V_gas_std;
    info.V_liq_std = V_liq_std;
    
    %% Step 4: Calculate Formation Volume Factors (Bo and Bg)
    
    try

        [~, Z_res_liq] = fugacitycoef_multicomp_liquid(comp, P_res, T_res,Pc, Tc, acentric, BIP, eos_type);
        [~, Z_res_vap] = fugacitycoef_multicomp(comp, P_res, T_res,Pc, Tc, acentric, BIP, eos_type);
        
        if Z_res_liq > 0 && Z_res_liq < Z_res_vap
            Z_res = Z_res_liq;
        else
            Z_res = Z_res_vap;
        end
        
        V_res_molar_EOS = Z_res * R * T_res / P_res;
        
        c_VT_res = calc_volume_translation(comp, P_res, T_res, Pc, Tc, acentric,vt_method, vt_params, n, eos_type, M_gmol,BIP);
        c_mix_res = sum(comp .* c_VT_res);
        V_res_molar = V_res_molar_EOS - c_mix_res;
        
        if V_res_molar <= 0
            V_res_molar = V_res_molar_EOS;
        end
        
    catch ME
        if verbose
            fprintf('EOS reservoir volume failed: %s\n', ME.message);
        end
        
        MW_res = sum(comp .* M_gmol);
        rho_res = 600;
        V_res_molar = MW_res / (rho_res * 1000);
    end
    
    if V_liq_std > 1e-15
        % Bo = V_res_molar / V_liq_molar;
        Bo = V_res_molar / ((1 - V_frac) * V_liq_molar);   % CORRECT
    else
        Bo = NaN;
    end
    
    V_std_gas = R * T_STO / P_STO;
    Bg = V_res_molar / V_std_gas;
    
    info.V_res_molar = V_res_molar;
    info.V_std_gas = V_std_gas;
    info.Bg = Bg;
    
    % Rs = GOR;
    

    
end


% Helper: Solve Rachford-Rice equation

function V = solve_RR(z, K)

    z = z(:);
    K = K(:);
    
    K_max = max(K);
    K_min = min(K);
    
    if K_max < 1
        V = 0;
        return;
    end
    if K_min > 1
        V = 1;
        return;
    end
    
    V_min = 1 / (1 - K_max);
    V_max = 1 / (1 - K_min);
    
    V_min = max(0, V_min + 1e-10);
    V_max = min(1, V_max - 1e-10);
    
    if V_min >= V_max
        V = 0.5;
        return;
    end
    
    tol = 1e-12;
    for iter = 1:100
        V = (V_min + V_max) / 2;
        f = sum(z .* (K - 1) ./ (1 + V * (K - 1)));
        
        if abs(f) < tol
            break;
        end
        
        if f > 0
            V_min = V;
        else
            V_max = V;
        end
    end
end


% Helper: Calculate volume translation

function c = calc_volume_translation(comp, press, temp, Pc, Tc, acentric, vt_method, vt_params, n, eos_type, M_gmol,BIP)

    c = zeros(n, 1);
    
    if isnumeric(vt_method) && length(vt_method) > 1
        c = vt_method(:) * 1e-6;
        if length(c) ~= n
            error('Volume translation vector length (%d) must match number of components (%d)', length(c), n);
        end
        return;
    end
    
    method = vt_method;
    if method == 0
        return;
    end
    
    switch method
        case 1  % Peneloux
            result = peneloux_volume_shift(Pc, Tc, acentric, comp, vt_params.Vc, eos_type);
            c = result.c_i * 1e-6;
            
        case 2  % Magoulas-Tassios
            [~, c, ~] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vt_params.Vc, vt_params.components);
            
        case 3  % Ungerer-Batut
            [~, c, ~] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vt_params.Vc);
            
        case 4  % Baled
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_params.components, comp, vt_params.Vc, eos_type);
            
        case 5  % Abudour
            if ~isfield(vt_params, 'Vc') || ~isfield(vt_params, 'components')
                error('Abudour method requires vt_params.Vc and vt_params.components');
            end
            Vc = vt_params.Vc;
            components = vt_params.components;
            % BIP_local = zeros(n);
            [~, c] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, M_gmol, components, BIP, eos_type);
      
        case 6  % Custom
            if ~isfield(vt_params, 'c_custom')
                error('Method 6 requires vt_params.c_custom');
            end
            c = vt_params.c_custom(:) * 1e-6;

        case 8  % Pina-Martinez (SRK)
            components = vt_params.components;
            result = pina_martinez_volume_shift(temp, Pc, Tc, acentric, components, comp);
            c = result.c_i * 1e-6;

        case 9  % Chen-Li (SRK)
            if isfield(vt_params, 'Zc') && ~isempty(vt_params.Zc)
                Zc = vt_params.Zc(:);
            else
                Zc = 0.2905 - 0.085 * acentric(:);
            end
            if isfield(vt_params, 'components') && ~isempty(vt_params.components)
                components = vt_params.components;
            else
                components = cell(n, 1);
            end
            [~, c, ~] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, Zc, components, false);

        case 10  % Baled (SRK)
            components = vt_params.components;
            Vc = vt_params.Vc;
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, components, comp, Vc, 'SRK');

        otherwise
            warning('Unknown VT method %d. Using no volume translation.', method);
    end

    c = c(:);    
 
end