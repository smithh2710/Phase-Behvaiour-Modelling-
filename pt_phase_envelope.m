function [env] = pt_phase_envelope(z, Pc, Tc, acentric, BIP, eos_type, comp_names, varargin)

p = inputParser;
addParameter(p, 'T_range', [], @isnumeric);
addParameter(p, 'tol', 1e-8, @isnumeric);
addParameter(p, 'maxiter', 25, @isnumeric);
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'T_res', NaN, @isnumeric);
addParameter(p, 'P_res', NaN, @isnumeric);
addParameter(p, 'M', [], @isnumeric);
addParameter(p, 'max_steps', 500, @isnumeric);
addParameter(p, 'T_crit', NaN, @isnumeric);
addParameter(p, 'P_crit', NaN, @isnumeric);
addParameter(p, 'ds_init', 0.12, @isnumeric);
addParameter(p, 'formulation', 'bubble', @ischar);
addParameter(p, 'beta', 0.5, @isnumeric);
addParameter(p, 'K_init', [], @isnumeric);
addParameter(p, 'T_init', NaN, @isnumeric);
addParameter(p, 'P_init', NaN, @isnumeric);
parse(p, varargin{:});
opts = p.Results;

z = z(:) / sum(z);
N = length(z);
is_dew = strcmpi(opts.formulation, 'dew');
is_quality = strcmpi(opts.formulation, 'quality');
beta = opts.beta;

if is_quality
    fprintf('Michelsen continuation (QUALITY formulation: V/F = %.2f)\n', beta);
elseif is_dew
    fprintf('Michelsen continuation (DEW formulation: y=z, sum(z/K)=1)\n');
else
    fprintf('Michelsen continuation (BUBBLE formulation: x=z, sum(Kz)=1)\n');
end

%% INITIALIZATION
has_K_init = ~isempty(opts.K_init) && ~isnan(opts.T_init) && ~isnan(opts.P_init);

if has_K_init
    K0 = opts.K_init(:);
    T0 = opts.T_init;
    P0 = opts.P_init;
    fprintf('  Using K-values from flash at T=%.1f C, P=%.1f bar\n', T0-273.15, P0/1e5);
    u0 = [log(K0); log(T0); log(P0)];

    [u0, conv0] = newton_solve(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
        N+1, log(T0), is_dew, is_quality, beta, opts.tol, opts.maxiter);
    if ~conv0
        [u0_try, conv0] = newton_solve([log(K0); log(T0); log(P0)], z, Pc, Tc, acentric, BIP, ...
            eos_type, comp_names, N+2, log(P0), is_dew, is_quality, beta, opts.tol, opts.maxiter);
        if conv0, u0 = u0_try; end
    end
else
    if isempty(opts.T_range)
        T_start = max(120, min(Tc) * 0.55);
    else
        T_start = opts.T_range(1);
    end

    K_w = wilson_K(T_start, 1e5, Pc, Tc, acentric);
    if is_dew
        P_w = 1e5 / max(sum(z ./ K_w), 1e-30);
    else
        P_w = sum(z .* K_w .* 1e5);
    end
    if P_w <= 0 || ~isfinite(P_w) || P_w > 1e9, P_w = 1e5; end
    u0 = [log(K_w); log(T_start); log(P_w)];

    [u0, conv0] = newton_solve(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
        N+1, log(T_start), is_dew, is_quality, beta, opts.tol, opts.maxiter);

    if ~conv0
        for T_s = [150, 180, 200, 250, max(120, min(Tc)*0.7)]
            K_w = wilson_K(T_s, 1e5, Pc, Tc, acentric);
            if is_dew
                P_w = 1e5 / max(sum(z ./ K_w), 1e-30);
            else
                P_w = sum(z .* K_w .* 1e5);
            end
            if P_w <= 0 || ~isfinite(P_w) || P_w > 1e9, continue; end
            u0 = [log(K_w); log(T_s); log(P_w)];
            [u0, conv0] = newton_solve(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
                N+1, log(T_s), is_dew, is_quality, beta, opts.tol, opts.maxiter);
            if conv0, break; end
        end
    end
end

if ~conv0
    warning('Could not initialize.'); env = empty_env(); return;
end
fprintf('  Start: T = %.1f C, P = %.1f bar\n', exp(u0(N+1))-273.15, exp(u0(N+2))/1e5);

%% CONTINUATION
if has_K_init
    fprintf('\n--- Forward continuation ---\n');
    [fwd_u, n_fwd] = run_continuation(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
        is_dew, is_quality, beta, +1, opts);

    fprintf('\n--- Backward T-stepping ---\n');
    [bwd_u, n_bwd] = run_tstep(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
        is_dew, is_quality, beta, -1, opts);

    all_u = [fliplr(bwd_u(:, 1:n_bwd)), u0, fwd_u(:, 1:n_fwd)];
    n_pts = n_bwd + 1 + n_fwd;
    fprintf('\n  Merged: %d bwd + 1 seed + %d fwd = %d total\n', n_bwd, n_fwd, n_pts);
else
    [fwd_u, n_fwd] = run_continuation(u0, z, Pc, Tc, acentric, BIP, eos_type, comp_names, ...
        is_dew, is_quality, beta, +1, opts);
    all_u = [u0, fwd_u(:, 1:n_fwd)];
    n_pts = 1 + n_fwd;
end

%% BUILD OUTPUT
T_all = exp(all_u(N+1, 1:n_pts))';
P_all = exp(all_u(N+2, 1:n_pts))';

T_crit = opts.T_crit;
P_crit = opts.P_crit;

if ~isnan(T_crit)
    i_bub = T_all <= T_crit;
    i_dew = T_all >= T_crit;
else
    [~, i_max] = max(P_all);
    i_bub = (1:n_pts)' <= i_max;
    i_dew = (1:n_pts)' >= i_max;
end

env.T_bub = T_all(i_bub); env.P_bub = P_all(i_bub);
env.T_dew = T_all(i_dew); env.P_dew = P_all(i_dew);
env.T_crit = T_crit; env.P_crit = P_crit;
env.T_all = T_all; env.P_all = P_all;

fprintf('\n=== Phase Envelope Summary ===\n');
fprintf('  Total: %d pts\n', n_pts);
if ~isnan(T_crit), fprintf('  Critical:  T = %.1f C, P = %.1f bar (input)\n', T_crit-273.15, P_crit/1e5); end

if opts.plot, do_plot(env, eos_type, opts); end
end


%% ================================================================
function [all_u, n_pts] = run_continuation(u0, z, Pc, Tc, acentric, BIP, ...
    eos_type, comp_names, is_dew, is_quality, beta, dir_sign, opts)

    N = length(z);
    max_steps = opts.max_steps;
    ds = opts.ds_init;
    ds_min = 0.0005;
    ds_max = 0.5;

    all_u = zeros(N+2, max_steps);
    n_pts = 0;
    u = u0;
    t_prev = [];
    n_fail = 0;

    for step = 1:max_steps
        J = numerical_jacobian(u, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
        [~, ~, V] = svd(J, 0);
        t = V(:, end);

        if ~isempty(t_prev)
            if dot(t, t_prev) < 0, t = -t; end
        else
            if dir_sign * t(N+1) < 0, t = -t; end
        end
        t_prev = t;

        [~, spec_idx] = max(abs(t));
        u_pred = u + ds * t;
        spec_val = u_pred(spec_idx);

        [u_new, conv, n_iter] = newton_solve(u_pred, z, Pc, Tc, acentric, BIP, ...
            eos_type, comp_names, spec_idx, spec_val, is_dew, is_quality, beta, opts.tol, opts.maxiter);

        if ~conv
            ds = ds * 0.35;
            if ds < ds_min
                n_fail = n_fail + 1;
                if n_fail > 10, fprintf('  Step %d: stalled\n', step); break; end
                ds = ds_min * 4;
            end
            continue;
        end
        n_fail = 0;

        n_pts = n_pts + 1;
        all_u(:, n_pts) = u_new;

        if mod(n_pts, 50) == 0 || n_pts <= 2
            fprintf('  Step %3d: T=%6.1f C, P=%7.1f bar\n', ...
                n_pts, exp(u_new(N+1))-273.15, exp(u_new(N+2))/1e5);
        end

        if n_iter <= 3, ds = min(ds * 1.4, ds_max);
        elseif n_iter >= 8, ds = max(ds * 0.5, ds_min); end

        u = u_new;

        if exp(u(N+2))/1e5 < 0.1 && n_pts > 5
            fprintf('  P < 0.1 bar, stopping.\n'); break;
        end
    end
end


%% ================================================================
function [all_u, n_pts] = run_tstep(u0, z, Pc, Tc, acentric, BIP, ...
    eos_type, comp_names, is_dew, is_quality, beta, dir_sign, opts)

    N = length(z);
    max_steps = opts.max_steps;
    dT_base = 3;

    all_u = zeros(N+2, max_steps);
    n_pts = 0;
    u_cur = u0;
    T_cur = exp(u0(N+1));
    n_consec_fail = 0;

    for step = 1:max_steps
        dT = dT_base;

        if ~isnan(opts.T_crit)
            dist = abs(T_cur - opts.T_crit);
            if dist < 20, dT = 1;
            elseif dist < 50, dT = 2; end
        end

        T_next = T_cur + dir_sign * dT;
        if T_next < 80 || T_next > 2000, break; end

        u_try = u_cur;
        u_try(N+1) = log(T_next);

        [u_new, conv, ~] = newton_solve(u_try, z, Pc, Tc, acentric, BIP, ...
            eos_type, comp_names, N+1, log(T_next), is_dew, is_quality, beta, opts.tol, opts.maxiter);

        if ~conv && dT > 0.5
            n_sub = 5; dT_sub = dT / n_sub;
            u_sub = u_cur; sub_ok = true;
            for ss = 1:n_sub
                T_sub = T_cur + dir_sign * ss * dT_sub;
                u_sub(N+1) = log(T_sub);
                [u_sub, cs, ~] = newton_solve(u_sub, z, Pc, Tc, acentric, BIP, ...
                    eos_type, comp_names, N+1, log(T_sub), is_dew, is_quality, beta, opts.tol, opts.maxiter);
                if ~cs, sub_ok = false; break; end
            end
            if sub_ok, conv = true; u_new = u_sub; end
        end

        if conv
            P_new = exp(u_new(N+2));
            if P_new/1e5 < 0.1 && n_pts > 3
                fprintf('  T-step: P < 0.1 bar at T=%.1f C, stopping.\n', T_next-273.15);
                break;
            end
            if P_new/1e5 > 2000
                fprintf('  T-step: P > 2000 bar at T=%.1f C, skipping.\n', T_next-273.15);
                n_consec_fail = n_consec_fail + 1;
                T_cur = T_next;
                if n_consec_fail > 5, break; end
                continue;
            end

            n_pts = n_pts + 1;
            all_u(:, n_pts) = u_new;
            u_cur = u_new;
            T_cur = T_next;
            n_consec_fail = 0;

            if mod(n_pts, 30) == 0
                fprintf('  T-step %3d: T=%6.1f C, P=%7.1f bar\n', n_pts, T_cur-273.15, P_new/1e5);
            end
        else
            n_consec_fail = n_consec_fail + 1;
            T_cur = T_next;
            if n_consec_fail > 5
                fprintf('  T-step: %d failures at T=%.1f C, stopping.\n', n_consec_fail, T_cur-273.15);
                break;
            end
        end
    end
    fprintf('  T-stepping: %d points\n', n_pts);
end


%% RESIDUAL
function f = eval_residual(u, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta)
    N = length(z);
    lnK = u(1:N); T = exp(u(N+1)); P = exp(u(N+2));
    K = exp(lnK);

    if is_quality
        denom = 1 + beta * (K - 1);
        denom = sign(denom) .* max(abs(denom), 1e-30);
        x = z ./ denom;
        y = K .* x;
        x = max(x, 1e-30); x = x / sum(x);
        y = max(y, 1e-30); y = y / sum(y);
    elseif is_dew
        y = z;
        S = sum(z ./ K);
        x = (z ./ K) / max(S, 1e-30);
        x = max(x, 1e-30); x = x / sum(x);
        y = max(y, 1e-30); y = y / sum(y);
    else
        x = z;
        S = sum(K .* z);
        y = K .* z / max(S, 1e-30);
        x = max(x, 1e-30); x = x / sum(x);
        y = max(y, 1e-30); y = y / sum(y);
    end

    try
        [phi_L, Z_L] = fugacitycoef_multicomp(x, P, T, Pc, Tc, acentric, BIP, eos_type, comp_names);
        [phi_V, Z_V] = fugacitycoef_multicomp(y, P, T, Pc, Tc, acentric, BIP, eos_type, comp_names);
        if numel(Z_L) > 1 && size(phi_L,2) > 1, [~,ii] = min(Z_L); phi_L = phi_L(:,ii); end
        if numel(Z_V) > 1 && size(phi_V,2) > 1, [~,ii] = max(Z_V); phi_V = phi_V(:,ii); end
    catch
        f = ones(N+1, 1) * 1e10; return;
    end

    phi_L = phi_L(:); phi_V = phi_V(:);
    f = zeros(N+1, 1);
    f(1:N) = lnK + log(max(phi_V, 1e-30)) - log(max(phi_L, 1e-30));

    if is_quality
        f(N+1) = sum(z .* (K - 1) ./ (1 + beta * (K - 1)));
    elseif is_dew
        f(N+1) = sum(z ./ K) - 1;
    else
        f(N+1) = sum(K .* z) - 1;
    end
end

%% NUMERICAL JACOBIAN
function J = numerical_jacobian(u, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta)
    N = length(z);
    eps_fd = 1e-6;
    f0 = eval_residual(u, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
    J = zeros(N+1, N+2);
    for j = 1:(N+2)
        u_p = u;
        h = eps_fd * max(abs(u(j)), 1e-3);
        u_p(j) = u(j) + h;
        fp = eval_residual(u_p, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
        J(:, j) = (fp - f0) / h;
    end
end

%% NEWTON SOLVER
function [u, converged, n_iter] = newton_solve(u, z, Pc, Tc, acentric, BIP, ...
    eos_type, comp_names, spec_idx, spec_val, is_dew, is_quality, beta, tol, maxiter)

    N = length(z);
    converged = false; n_iter = maxiter;
    eps_fd = 1e-6;

    for iter = 1:maxiter
        f_red = eval_residual(u, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
        f_full = [f_red; u(spec_idx) - spec_val];

        if max(abs(f_full)) < tol
            converged = true; n_iter = iter; return;
        end
        if any(~isfinite(f_full)), return; end

        J_red = zeros(N+1, N+2);
        for j = 1:(N+2)
            u_p = u;
            h = eps_fd * max(abs(u(j)), 1e-3);
            u_p(j) = u(j) + h;
            fp = eval_residual(u_p, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
            J_red(:, j) = (fp - f_red) / h;
        end

        J_spec = zeros(1, N+2); J_spec(spec_idx) = 1.0;
        J_full = [J_red; J_spec];

        rc = rcond(J_full);
        if rc < 1e-18 || isnan(rc), J_full = J_full + 1e-10 * eye(N+2); end

        du = -J_full \ f_full;
        if max(abs(du)) > 2.0, du = du * (2.0 / max(abs(du))); end

        alpha = 1.0;
        f_norm_old = norm(f_full);
        for back = 1:10
            u_try = u + alpha * du;
            if exp(u_try(N+1)) < 50 || exp(u_try(N+1)) > 3000 || exp(u_try(N+2)) < 1 || exp(u_try(N+2)) > 1e9
                alpha = alpha * 0.5; continue;
            end
            f_try = eval_residual(u_try, z, Pc, Tc, acentric, BIP, eos_type, comp_names, is_dew, is_quality, beta);
            f_try_full = [f_try; u_try(spec_idx) - spec_val];
            if norm(f_try_full) < f_norm_old || alpha < 0.02, break; end
            alpha = alpha * 0.5;
        end
        u = u + alpha * du;
    end
end

function K = wilson_K(T, P, Pc, Tc, acentric)
    K = (Pc / P) .* exp(5.373 * (1 + acentric) .* (1 - Tc / T));
    K = max(K, 1e-30); K = min(K, 1e30);
end

function env = empty_env()
    env.T_bub = []; env.P_bub = [];
    env.T_dew = []; env.P_dew = [];
    env.T_all = []; env.P_all = [];
    env.T_crit = NaN; env.P_crit = NaN;
end

function do_plot(env, eos_type, opts)
    figure('Position', [100 100 700 550]);
    hold on;
    if ~isempty(env.T_bub)
        plot(env.T_bub-273.15, env.P_bub/1e5, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Bubble');
    end
    if ~isempty(env.T_dew)
        plot(env.T_dew-273.15, env.P_dew/1e5, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Dew');
    end
    if ~isnan(env.T_crit)
        plot(env.T_crit-273.15, env.P_crit/1e5, 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'm', ...
            'DisplayName', sprintf('Critical (%.0f°C, %.0f bar)', env.T_crit-273.15, env.P_crit/1e5));
    end
    if ~isnan(opts.T_res) && ~isnan(opts.P_res)
        plot(opts.T_res-273.15, opts.P_res/1e5, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
            'LineWidth', 1.5, 'DisplayName', sprintf('Reservoir (%.0f°C, %.0f bar)', opts.T_res-273.15, opts.P_res/1e5));
    end
    hold off;
    xlabel('Temperature (°C)'); ylabel('Pressure (bar)');
    title(sprintf('PT Phase Envelope — %s EOS', eos_type));
    legend('Location', 'best', 'FontSize', 8); grid on; set(gca, 'FontSize', 11);
end