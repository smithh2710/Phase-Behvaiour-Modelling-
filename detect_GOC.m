function [GOC_depth, GOC_pressure, GOC_temp, GOC_comp, GOC_type, info] = detect_GOC(...
    comp_ref, press_ref, temp_ref, h_ref, depth_range, ...
    Pc, Tc, acentric, BIP, M_gmol, ...
    dTdh, Cp_coeffs, H_ig_ref, vt_method, vt_params, varargin)

    p = inputParser;
    addParameter(p, 'scan_step', 5, @isnumeric);
    addParameter(p, 'tol_depth', 0.5, @isnumeric);
    addParameter(p, 'tol_DeltaK', 1e-6, @isnumeric);
    addParameter(p, 'max_refine', 40, @isnumeric);
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'grading_model', 'hasse', @ischar);
    addParameter(p, 'tau', 4.0, @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;
    tau = opts.tau;

    comp_ref = comp_ref(:);
    n = length(comp_ref);
    tol_sat = 1e-8;
    maxiter_sat = 500;

    h_min = min(depth_range);
    h_max = max(depth_range);

    if isnumeric(vt_method) && isscalar(vt_method)
        if (vt_method >= 7 && vt_method <= 10) || vt_method == 12
            eos_type = 'SRK';
        else
            eos_type = 'PR';
        end
    else
        eos_type = 'PR';
    end

    if isfield(vt_params, 'components') && ~isempty(vt_params.components)
        components = vt_params.components;
    else
        components = cell(n, 1);
    end

    GOC_depth = NaN;
    GOC_pressure = NaN;
    GOC_temp = NaN;
    GOC_comp = comp_ref;
    GOC_type = 'none';

    function [DK, Dp, P_h, T_h, z_h, P_sat, sat_type, K_vals] = eval_depth(h)
        switch lower(opts.grading_model)
            case {'hasse', 'isothermal'}
                dTdh_use = dTdh;
                if strcmpi(opts.grading_model, 'isothermal'), dTdh_use = 0; end
                [z_h, P_h, T_h, ~, ~] = main_hasse(h, h_ref, comp_ref, press_ref, temp_ref, ...
                    dTdh_use, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
            case 'kempers'
                [z_h, P_h, T_h, ~, ~] = main_kempers(h, h_ref, comp_ref, press_ref, temp_ref, ...
                    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
            case 'firoozabadi'
                [z_h, P_h, T_h, ~, ~] = main_firoozabadi(h, h_ref, comp_ref, press_ref, temp_ref, ...
                    dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
        end
        z_h = z_h(:);

        [P_bub_calc, K_bub, DK_bub] = solve_bubble(z_h, P_h, T_h);
        [P_dew_calc, K_dew, DK_dew] = solve_dew(z_h, P_h, T_h);

        if ~isnan(P_bub_calc) && ~isnan(P_dew_calc)
            if P_bub_calc >= P_dew_calc
                sat_type = 'Bubble'; P_sat = P_bub_calc; K_vals = K_bub; DK = DK_bub;
            else
                sat_type = 'Dew'; P_sat = P_dew_calc; K_vals = K_dew; DK = DK_dew;
            end
        elseif ~isnan(P_bub_calc)
            sat_type = 'Bubble'; P_sat = P_bub_calc; K_vals = K_bub; DK = DK_bub;
        elseif ~isnan(P_dew_calc)
            sat_type = 'Dew'; P_sat = P_dew_calc; K_vals = K_dew; DK = DK_dew;
        else
            sat_type = 'Unknown'; P_sat = NaN; K_vals = ones(n,1); DK = NaN;
        end

        Dp = (P_h - P_sat) / P_h;
    end

    function [P_calc, K, DK] = solve_bubble(z, P, T)
        P_calc = NaN; K = ones(n,1); DK = NaN;
        try
            [P_calc, comp_vap] = pressbub_multicomp_newton(z, P, T, Pc, Tc, acentric, BIP, tol_sat, maxiter_sat, eos_type, components);
            if ~isreal(P_calc) || P_calc <= 0 || ~isfinite(P_calc)
                error('invalid');
            end
            K = comp_vap(:) ./ z;
            K(z < 1e-15) = 1;
            if max(abs(K - 1)) < 1e-3
                P_calc = NaN; K = ones(n,1); DK = NaN; return;
            end
            DK = sum(log(K).^2);
        catch
            try
                [P_calc, comp_vap] = pressbub_multicomp_ss(z, P, T, Pc, Tc, acentric, BIP, tol_sat, maxiter_sat, eos_type, components);
                if ~isreal(P_calc) || P_calc <= 0 || ~isfinite(P_calc)
                    P_calc = NaN; K = ones(n,1); DK = NaN; return;
                end
                K = comp_vap(:) ./ z;
                K(z < 1e-15) = 1;
                if max(abs(K - 1)) < 1e-3
                    P_calc = NaN; K = ones(n,1); DK = NaN; return;
                end
                DK = sum(log(K).^2);
            catch
                P_calc = NaN; K = ones(n,1); DK = NaN;
            end
        end
    end

    function [P_calc, K, DK] = solve_dew(z, P, T)
        P_calc = NaN; K = ones(n,1); DK = NaN;
        try
            [P_calc, comp_liq] = pressdew_multicomp_newton(z, P, T, Pc, Tc, acentric, BIP, tol_sat, maxiter_sat, eos_type, components);
            if ~isreal(P_calc) || P_calc <= 0 || ~isfinite(P_calc)
                error('invalid');
            end
            K = z ./ comp_liq(:);
            K(z < 1e-15) = 1;
            K(comp_liq < 1e-15) = 1;
            if max(abs(K - 1)) < 1e-3
                P_calc = NaN; K = ones(n,1); DK = NaN; return;
            end
            DK = sum(log(K).^2);
        catch
            try
                [P_calc, comp_liq] = pressdew_multicomp_ss(z, P, T, Pc, Tc, acentric, BIP, tol_sat, maxiter_sat, eos_type, components);
                if ~isreal(P_calc) || P_calc <= 0 || ~isfinite(P_calc)
                    P_calc = NaN; K = ones(n,1); DK = NaN; return;
                end
                K = z ./ comp_liq(:);
                K(z < 1e-15) = 1;
                K(comp_liq < 1e-15) = 1;
                if max(abs(K - 1)) < 1e-3
                    P_calc = NaN; K = ones(n,1); DK = NaN; return;
                end
                DK = sum(log(K).^2);
            catch
                P_calc = NaN; K = ones(n,1); DK = NaN;
            end
        end
    end

    %% Phase 1: Coarse scan
    scan_depths = h_min : opts.scan_step : h_max;
    if scan_depths(end) < h_max
        scan_depths = [scan_depths, h_max];
    end
    n_scan = length(scan_depths);

    scan_DK   = NaN(n_scan, 1);
    scan_Dp   = NaN(n_scan, 1);
    scan_P    = NaN(n_scan, 1);
    scan_T    = NaN(n_scan, 1);
    scan_Psat = NaN(n_scan, 1);
    scan_type = cell(n_scan, 1);
    scan_comp = zeros(n_scan, n);

    if opts.verbose
        fprintf('=== GOC Detection: %s, EOS=%s, VT=%d ===\n', upper(opts.grading_model), eos_type, vt_method);
        fprintf('%-8s %-10s %-10s %-10s %-10s %-10s\n', 'Depth', 'P_res', 'P_sat', 'Delta_p', 'Delta_K', 'Type');
        fprintf('%s\n', repmat('-', 1, 62));
    end

    sat_bracket = [];
    unsat_bracket = [];

    for i = 1:n_scan
        h = scan_depths(i);
        try
            [DK, Dp, P_h, T_h, z_h, P_sat, sat_type, ~] = eval_depth(h);
        catch
            continue;
        end

        scan_DK(i) = DK;
        scan_Dp(i) = Dp;
        scan_P(i) = P_h;
        scan_T(i) = T_h;
        scan_Psat(i) = P_sat;
        scan_type{i} = sat_type;
        scan_comp(i,:) = z_h';

        if opts.verbose
            fprintf('%-8.1f %-10.2f %-10.2f %-10.4f %-10.4f %-10s\n', ...
                h, P_h/1e5, P_sat/1e5, Dp, DK, sat_type);
        end

        if i < 2, continue; end
        if isnan(scan_DK(i)) || isnan(scan_Dp(i)), continue; end
        if isnan(scan_DK(i-1)) || isnan(scan_Dp(i-1)), continue; end

        if isempty(sat_bracket) && scan_Dp(i) * scan_Dp(i-1) < 0
            sat_bracket = [scan_depths(i-1), scan_depths(i)];
            if opts.verbose
                fprintf('  >> Saturated GOC bracket: [%.1f, %.1f] m\n', sat_bracket(1), sat_bracket(2));
            end
        end

        if isempty(unsat_bracket)
            if DK < opts.tol_DeltaK
                unsat_bracket = [scan_depths(i-1), scan_depths(i)];
                if opts.verbose
                    fprintf('  >> Undersaturated GOC bracket: [%.1f, %.1f] m (DK=%.2e)\n', ...
                        unsat_bracket(1), unsat_bracket(2), DK);
                end
            end

            if ~isempty(scan_type{i-1}) && ~strcmp(sat_type, scan_type{i-1})
                if i > 2 && ~isempty(scan_type{i-2}) && strcmp(scan_type{i-2}, scan_type{i-1})
                    unsat_bracket = [scan_depths(i-1), scan_depths(i)];
                    if opts.verbose
                        fprintf('  >> Undersaturated GOC bracket: [%.1f, %.1f] m (%s -> %s)\n', ...
                            unsat_bracket(1), unsat_bracket(2), scan_type{i-1}, sat_type);
                    end
                end
            end

            if i > 2 && ~isnan(scan_DK(i-2))
                if scan_DK(i) > scan_DK(i-1) && scan_DK(i-1) < scan_DK(i-2)
                    unsat_bracket = [scan_depths(i-2), scan_depths(i)];
                    if opts.verbose
                        fprintf('  >> Undersaturated GOC bracket: [%.1f, %.1f] m (DK minimum)\n', ...
                            unsat_bracket(1), unsat_bracket(2));
                    end
                end
            end
        end

        if ~isempty(sat_bracket) && ~isempty(unsat_bracket)
            break;
        end
    end

    %% Store scan info
    info.scan_depths = scan_depths(:);
    info.scan_DK = scan_DK;
    info.scan_Dp = scan_Dp;
    info.scan_P = scan_P;
    info.scan_T = scan_T;
    info.scan_Psat = scan_Psat;
    info.scan_type = scan_type;

    %% Phase 2: Refine
    sat_GOC_h = NaN;  sat_GOC_P = NaN; sat_GOC_T = NaN; sat_GOC_z = comp_ref;
    unsat_GOC_h = NaN; unsat_GOC_P = NaN; unsat_GOC_T = NaN; unsat_GOC_z = comp_ref;

    if ~isempty(sat_bracket)
        h_lo = sat_bracket(1); h_hi = sat_bracket(2);
        [~, Dp_lo, ~, ~, ~, ~, ~, ~] = eval_depth(h_lo);
        for iter = 1:opts.max_refine
            h_mid = (h_lo + h_hi) / 2;
            [~, Dp_mid, P_mid, T_mid, z_mid, ~, ~, ~] = eval_depth(h_mid);
            if abs(h_hi - h_lo) < opts.tol_depth || abs(Dp_mid) < 1e-6
                break;
            end
            if Dp_lo * Dp_mid < 0
                h_hi = h_mid;
            else
                h_lo = h_mid; Dp_lo = Dp_mid;
            end
        end
        sat_GOC_h = h_mid; sat_GOC_P = P_mid; sat_GOC_T = T_mid; sat_GOC_z = z_mid;
    end

    if ~isempty(unsat_bracket)
        h_lo = unsat_bracket(1); h_hi = unsat_bracket(2);
        [DK_lo, ~, ~, ~, ~, ~, type_lo, ~] = eval_depth(h_lo);
        [DK_hi, ~, ~, ~, ~, ~, type_hi, ~] = eval_depth(h_hi);

        for iter = 1:opts.max_refine
            h_mid = (h_lo + h_hi) / 2;
            [DK_mid, ~, P_mid, T_mid, z_mid, ~, type_mid, ~] = eval_depth(h_mid);
            if abs(h_hi - h_lo) < opts.tol_depth || DK_mid < opts.tol_DeltaK
                break;
            end
            if ~strcmp(type_lo, type_mid)
                h_hi = h_mid; DK_hi = DK_mid; type_hi = type_mid;
            elseif ~strcmp(type_hi, type_mid)
                h_lo = h_mid; DK_lo = DK_mid; type_lo = type_mid;
            else
                if DK_lo < DK_hi
                    h_hi = h_mid; DK_hi = DK_mid;
                else
                    h_lo = h_mid; DK_lo = DK_mid;
                end
            end
        end
        unsat_GOC_h = h_mid; unsat_GOC_P = P_mid; unsat_GOC_T = T_mid; unsat_GOC_z = z_mid;
    end

    %% Phase 3: Pick the deeper GOC
    if ~isnan(sat_GOC_h) && ~isnan(unsat_GOC_h)
        if sat_GOC_h >= unsat_GOC_h
            GOC_depth = sat_GOC_h; GOC_pressure = sat_GOC_P;
            GOC_temp = sat_GOC_T; GOC_comp = sat_GOC_z; GOC_type = 'saturated';
        else
            GOC_depth = unsat_GOC_h; GOC_pressure = unsat_GOC_P;
            GOC_temp = unsat_GOC_T; GOC_comp = unsat_GOC_z; GOC_type = 'undersaturated';
        end
    elseif ~isnan(sat_GOC_h)
        GOC_depth = sat_GOC_h; GOC_pressure = sat_GOC_P;
        GOC_temp = sat_GOC_T; GOC_comp = sat_GOC_z; GOC_type = 'saturated';
    elseif ~isnan(unsat_GOC_h)
        GOC_depth = unsat_GOC_h; GOC_pressure = unsat_GOC_P;
        GOC_temp = unsat_GOC_T; GOC_comp = unsat_GOC_z; GOC_type = 'undersaturated';
    end

    info.sat_GOC_h = sat_GOC_h;
    info.unsat_GOC_h = unsat_GOC_h;
    info.GOC_depth = GOC_depth;
    info.GOC_type = GOC_type;

    if opts.verbose
        fprintf('\n--- Result ---\n');
        if ~isnan(sat_GOC_h), fprintf('  Saturated GOC:      %.2f m\n', sat_GOC_h); end
        if ~isnan(unsat_GOC_h), fprintf('  Undersaturated GOC: %.2f m\n', unsat_GOC_h); end
        if strcmp(GOC_type, 'none')
            fprintf('  No GOC found in [%.0f, %.0f] m\n', h_min, h_max);
        else
            fprintf('  Selected (%s): h=%.2f m, P=%.2f bar, T=%.1f K\n', ...
                GOC_type, GOC_depth, GOC_pressure/1e5, GOC_temp);
        end
    end
end