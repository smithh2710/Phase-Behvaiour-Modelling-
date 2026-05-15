function [rho, V_m, Z, c_mix] = calculate_density(comp, press, temp, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params)

    R = 8.3144598;
    
    if nargin < 10
        vt_params = struct();
    end
    if nargin < 9
        vt_method = 0;
    end
    
    comp = comp(:);
    comp = comp / sum(comp);
    n = length(comp); 
    
    [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol);
    
    % P-dependent methods: solve EOS internally, use struct directly
    if vt_type == 9
        [~, ~, mix_cl] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Zc, vt_opts.components, BIP, false);
        c_mix = mix_cl.c_mix;
        V_m_eos =  mix_cl.V_EOS;
        Z = mix_cl.Z;

    elseif vt_type == 5
        [~, ~, mix_abu] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Vc, M_gmol, vt_opts.components, BIP, 'PR');
        c_mix = mix_abu.c_mix;
        V_m_eos = mix_abu.V_EOS;
        Z = mix_abu.Z;

    else
        % All other methods: solve EOS here, get c from VT functions
        [~, Z] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP, eos_type);
        V_m_eos = Z * R * temp / press;
        
        switch vt_type
            case {0, 7}
                c_mix = 0;

            % Struct-based methods (mixture c_mix directly, all in m³/mol)
            case 2
                [~, ~, mix] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vt_opts.Vc, vt_opts.components);
                c_mix = mix.c_mix;

            case 3
                [~, ~, mix] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vt_opts.Vc);
                c_mix = - mix.c_mix;

            case {4, 10}
                if vt_type == 4
                    eos_baled = 'PR';
                else
                    eos_baled = 'SRK';
                end
                [~, ~, mix] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_opts.components, comp, vt_opts.Vc, eos_baled);
                c_mix = mix.c_mix;

            % Component c_i methods (all return c in m³/mol now)
            otherwise
                c_i = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
                c_mix = sum(comp .* c_i);
        end
    end
    
    V_m = V_m_eos - c_mix;
    M_mix = sum(comp .* M_gmol) / 1000;
    rho = M_mix / V_m;

end


function c = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts)

    n = length(comp);

    switch vt_type
        case 1
            [~, c] = peneloux_volume_shift(Pc, Tc, acentric, comp, vt_opts.Vc, 'PR');
        case 8
            [~, c] = pina_martinez_volume_shift(temp, Pc, Tc, acentric, vt_opts.components, comp);
        case 11
            [~, c] = jhaveri_youngren_volume_shift(Pc / 1e5, Tc, acentric, M_gmol);
            c = c * 1e-6;
        case {6, 12}
            c = vt_opts.c_direct;
        otherwise
            c = zeros(n, 1);
    end
    
    c = c(:);
end


function [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol)

    R = 8.3144598;
    vt_opts = struct();
    
    if isfield(vt_params, 'Vc') && ~isempty(vt_params.Vc)
        vt_opts.Vc = vt_params.Vc(:);
    else
        Zc_est = 0.2905 - 0.085 * acentric;
        vt_opts.Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
    end
    
    if isfield(vt_params, 'Zc') && ~isempty(vt_params.Zc)
        vt_opts.Zc = vt_params.Zc(:);
    else
        vt_opts.Zc = 0.2905 - 0.085 * acentric(:);
    end
    
    if isfield(vt_params, 'components') && ~isempty(vt_params.components)
        vt_opts.components = vt_params.components;
    else
        vt_opts.components = cell(n, 1);
    end
    
    vt_opts.Pc = Pc;
    vt_opts.Tc = Tc;
    vt_opts.acentric = acentric;
    vt_opts.M_gmol = M_gmol;
    vt_opts.temp = temp;
    vt_opts.n = n;
    
    if isempty(vt_method)
        vt_type = 0;
        eos_type = 'PR';
        return;
    end
    
    if isnumeric(vt_method) && length(vt_method) == n && n > 1
        vt_type = -1;
        eos_type = 'PR';
        vt_opts.c_direct = vt_method(:) * 1e-6;
        return;
    end
    
    if ischar(vt_method) || isstring(vt_method)
        vt_type = string_to_method_number(vt_method);
    else
        vt_type = vt_method;
    end
    
    if vt_type >= 0 && vt_type <= 6 || vt_type == 11
        eos_type = 'PR';
    elseif vt_type >= 7 && vt_type <= 10 || vt_type == 12
        eos_type = 'SRK';
    else
        vt_type = 0;
        eos_type = 'PR';
    end
    
    if vt_type == 6 || vt_type == 12
        if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
            vt_opts.c_direct = vt_params.c_custom(:) * 1e-6;
        else
            error('For vt_method = %d, provide vt_params.c_custom [cm³/mol]', vt_type);
        end
    end

end


function num = string_to_method_number(name)
    name = lower(char(name));
    switch name
        case {'pr', '0'}, num = 0;
        case {'peneloux', '1'}, num = 1;
        case {'magoulas', 'mt', '2'}, num = 2;
        case {'ungerer', 'ub', '3'}, num = 3;
        case {'baled_pr', '4'}, num = 4;
        case {'abudour', '5'}, num = 5;
        case {'custom', '6'}, num = 6;
        case {'jhaveri', 'jy', '11'}, num = 11;
        case {'srk', '7'}, num = 7;
        case {'pina', 'pm', '8'}, num = 8;
        case {'chen_li', 'cl', '9'}, num = 9;
        case {'baled_srk', '10'}, num = 10;   
        case {'srk_custom', 'srk_peneloux', '12'}, num = 12;
        otherwise, num = 0;
    end
end