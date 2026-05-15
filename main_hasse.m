function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_hasse(h_target, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)

    if nargin < 15
        vt_params = struct();
    end
    if nargin < 14
        vt_method = 0;
    end

    n = length(comp_ref);
    step_size = 1;
    
    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    H_ig_ref = H_ig_ref(:);

    [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, T_ref, Pc, Tc, acentric, M_gmol);
    vt_opts.BIP = BIP;

    tol = 1e-10;

    components = vt_opts.components;

    
    maxiter = 1500;

    delta_h_total = h_target - h_ref;
    n_steps = ceil(abs(delta_h_total) / step_size);
    
    if n_steps == 0
        n_steps = 1;
    end
    
    dh = delta_h_total / n_steps;

    h_current = h_ref;
    P_current = P_ref;
    T_current = T_ref;
    comp_current = comp_ref;

    for step = 1:n_steps
        h_next = h_current + dh;
        
        [comp_next, P_next, T_next] = hasse_single_step(...
            h_next, h_current, comp_current, P_current, T_current, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts, components);
        
        h_current = h_next;
        P_current = P_next;
        T_current = T_next;
        comp_current = comp_next;
    end

    comp_h = comp_current;
    press_h = P_current;
    temp_h = T_current;

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

    try 
        pressbub_ini_h = press_h;
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, pressbub_ini_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
        if ~isreal(pressbub_h) || pressbub_h <= 0 || ~isfinite(pressbub_h)
           pressbub_h = NaN;
        end
    catch 
       pressbub_h = NaN; 
    end
    if isnan(pressbub_h) || ~isreal(pressbub_h) || pressbub_h <= 0 || pressbub_h > press_h
        try
            [pressbub_h, ~] = pressbub_multicomp_ss(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
        catch
            pressbub_h = NaN;
        end
    end
    if ~isnan(pressbub_h) && (pressbub_h <= 0 || pressbub_h > press_h)
        pressbub_h = NaN;
    end

    try 
        pressdew_ini_h = 250e5;
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
        if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
           pressdew_h = NaN;
        end
    catch 
       pressdew_h = NaN; 
    end
    if isnan(pressdew_h) || ~isreal(pressdew_h) || pressdew_h <= 0 || pressdew_h > press_h
        try
            [pressdew_h, ~] = pressdew_multicomp_ss(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type, components);
        catch
            pressdew_h = NaN;
        end
    end
    if ~isnan(pressdew_h) && (pressdew_h <= 0 || pressdew_h > press_h)
        pressdew_h = NaN;
    end

end


function [comp_h, press_h, temp_h] = hasse_single_step(h, h_ref, comp_ref, P, T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts, components)

    R = 8.3144598;
    g = 9.80665;
    n = length(comp_ref);

    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);

    delta_h = h - h_ref;
    delta_T = dTdh * delta_h;
    temp_h = T + delta_T;

    [c_ref, fugcoef_internal] = calc_vt(comp_ref, P, T, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);

    if ~isempty(fugcoef_internal)
        fugcoef_ref_eos = fugcoef_internal;
    else
        [fugcoef_ref_eos, ~] = fugacitycoef_multicomp(comp_ref, P, T, Pc, Tc, acentric, BIP, eos_type, components);
    end
    
    if vt_type ~= 0 && vt_type ~= 7
        fugcoef_ref = fugcoef_ref_eos .* exp(-c_ref * P / (R * T));
    else
        fugcoef_ref = fugcoef_ref_eos;
    end
    
    f_ref = fugcoef_ref .* comp_ref * P;

    [H_mix_specific, H_partial_specific, ~, ~, ~] = calculate_absolute_enthalpy(T, P, comp_ref, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type, components);

    term2 = zeros(n, 1);
    term3 = zeros(n, 1);

    for i = 1:n
        M_i_kg = M_gmol(i) / 1000;
        term2(i) = M_i_kg * g * delta_h / (R * T);
        
        enthalpy_diff = H_mix_specific - H_partial_specific(i);
        term3(i) = M_gmol(i) * enthalpy_diff * delta_T / (R * T^2);
    end

    params = struct();
    params.f_ref = f_ref;
    params.term2 = term2;
    params.term3 = term3;
    params.T_ref = T;
    params.Pc = Pc;
    params.Tc = Tc;
    params.acentric = acentric;
    params.BIP = BIP;
    params.M_gmol = M_gmol;
    params.n = n;
    params.vt_type = vt_type;
    params.eos_type = eos_type;
    params.vt_opts = vt_opts;
    params.R = R;
    params.components = components;

    initial_guess = [comp_ref; P];

    options = optimoptions('fsolve', ...
        'Display', 'none', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxIterations', 200);

    [solution, ~, ~] = fsolve(@(x) residual_hasse(x, params), initial_guess, options);

    comp_h = solution(1:n);
    press_h = solution(end);

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

end


function F = residual_hasse(x, params)

    f_ref = params.f_ref;
    term2 = params.term2;
    term3 = params.term3;
    T_ref = params.T_ref;
    Pc = params.Pc;
    Tc = params.Tc;
    acentric = params.acentric;
    BIP = params.BIP;
    M_gmol = params.M_gmol;
    n = params.n;
    vt_type = params.vt_type;
    eos_type = params.eos_type;
    vt_opts = params.vt_opts;
    R = params.R;
    components = params.components;

    z_h = x(1:n);
    P_h = x(end);

    z_h = max(z_h, 1e-15);
    z_h_norm = z_h / sum(z_h);
    P_h = max(P_h, 1e2);

    if vt_type ~= 0 && vt_type ~= 7
        [c_h, fugcoef_internal] = calc_vt(z_h_norm, P_h, T_ref, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
        if ~isempty(fugcoef_internal)
            fugcoef_h_eos = fugcoef_internal;
        else
            [fugcoef_h_eos, ~] = fugacitycoef_multicomp(z_h_norm, P_h, T_ref, Pc, Tc, acentric, BIP, eos_type, components);
        end
        fugcoef_h = fugcoef_h_eos .* exp(-c_h * P_h / (R * T_ref));
    else
        [fugcoef_h_eos, ~] = fugacitycoef_multicomp(z_h_norm, P_h, T_ref, Pc, Tc, acentric, BIP, eos_type, components);
        fugcoef_h = fugcoef_h_eos;
    end

    F = zeros(n+1, 1);

    for i = 1:n
        f_i_h = fugcoef_h(i) * z_h_norm(i) * P_h;
        F(i) = log(f_i_h / f_ref(i)) - term2(i) + term3(i);
    end

    F(n+1) = sum(z_h) - 1;

end


function [c, fugcoef] = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts)

    n = length(comp);
    comp = comp(:);
    fugcoef = [];
    
    switch vt_type
        
        case 0
            c = zeros(n, 1);
            
        case 1
            result = peneloux_volume_shift(Pc, Tc, acentric, comp, vt_opts.Vc, 'PR');
            c = result.c_i * 1e-6;
            
        case 2
            [~, c, ~] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vt_opts.Vc, vt_opts.components);
            
        case 3
            [~, c, ~] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vt_opts.Vc);
            
        case 4
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_opts.components, comp, vt_opts.Vc, 'PR');
            
        case 5
            [~, c, mix_abu] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Vc, M_gmol, vt_opts.components, vt_opts.BIP, 'PR');
            fugcoef = mix_abu.fugcoef;
            
        case 6
            c = vt_opts.c_direct;
            
        case 7
            c = zeros(n, 1);
            
        case 8
            result = pina_martinez_volume_shift(temp, Pc, Tc, acentric, vt_opts.components, comp);
            c = result.c_i * 1e-6;
            
        case 9
            [~, c, mix_cl] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Zc, vt_opts.components, vt_opts.BIP, false);
            fugcoef = - mix_cl.fugcoef;
            c = -c ; 
            
        case 10
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_opts.components, comp, vt_opts.Vc, 'SRK');
            
        case -1
            c = vt_opts.c_direct;

        case 11
            [~, c_cm3] = jhaveri_youngren_volume_shift(Pc / 1e5, Tc, acentric, M_gmol);
            c = c_cm3 * 1e-6;

        case 12
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
        warning('Unknown VT method %d. Defaulting to PR with no VT.', vt_type);
        vt_type = 0;
        eos_type = 'PR';
    end
    
    if vt_type == 6
        if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
            c_custom = vt_params.c_custom(:);
            if length(c_custom) ~= n
                error('c_custom must have %d elements (one per component)', n);
            end
            vt_opts.c_direct = c_custom * 1e-6;
        else
            error('For vt_method = 6, you must provide vt_params.c_custom [cm³/mol]');
        end
    end

    if vt_type == 12
        if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
            c_custom = vt_params.c_custom(:);
            if length(c_custom) ~= n
                error('c_custom must have %d elements (one per component)', n);
            end
            vt_opts.c_direct = c_custom * 1e-6;
        else
            error('For vt_method = 12, you must provide vt_params.c_custom [cm³/mol]');
        end
    end
    
end


function num = string_to_method_number(name)

    name = lower(char(name));
    
    switch name
        case {'pr', 'pr_novt', '0'}
            num = 0;
        case {'pr_peneloux', 'peneloux', 'pen', '1'}
            num = 1;
        case {'pr_mt', 'magoulas_tassios', 'magoulas', 'mt', '2'}
            num = 2;
        case {'pr_ub', 'ungerer_batut', 'ungerer', 'ub', '3'}
            num = 3;
        case {'pr_baled', 'baled_pr', '4'}
            num = 4;
        case {'pr_abudour', 'abudour', 'abu', '5'}
            num = 5;
        case {'pr_custom', 'custom', '6'}
            num = 6;
        case {'pr_jy', 'jhaveri_youngren', 'jhaveri', 'jy', '11'}
            num = 11;
        case {'srk', 'srk_novt', '7'}
            num = 7;
        case {'srk_pm', 'pina_martinez', 'pina', 'pm', '8'}
            num = 8;
        case {'srk_cl', 'chen_li', 'chen', 'cl', '9'}
            num = 9;
        case {'srk_baled', 'baled_srk', '10'}
            num = 10;
        case {'srk_custom', 'srk_peneloux', 'srk_pen', '12'}
            num = 12;
        otherwise
            warning('Unknown VT method: %s. Defaulting to PR with no VT.', name);
            num = 0;
    end
end