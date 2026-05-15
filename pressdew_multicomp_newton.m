function [pressd, comp_liq, K_out] = pressdew_multicomp_newton(comp_vap, pressd_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type, components, K_init)

if nargin < 10 || isempty(eos_type), eos_type = 'PR'; end
if nargin < 11, components = {}; end

ncomp = length(comp_vap);
comp_vap = comp_vap(:);

have_K_init = (nargin >= 12) && ~isempty(K_init) && all(K_init > 0) && all(isfinite(K_init));

if have_K_init
    K = K_init(:);
    P = pressd_ini;
else
    K = wilsoneq(pressd_ini, temp, pressc, tempc, acentric);
    P = pressd_ini;

    n_ss = min(15, round(maxiter * 0.15));
    for i = 1:n_ss
        x = comp_vap ./ K;
        sx = sum(x);
        x = max(x, 1e-20);
        x = x / sum(x);

        [phi_V, ~] = fugacitycoef_multicomp_vapor(comp_vap, P, temp, pressc, tempc, acentric, BIP, eos_type, components);
        [phi_L, ~] = fugacitycoef_multicomp_liquid(x, P, temp, pressc, tempc, acentric, BIP, eos_type, components);

        K = phi_L(:) ./ phi_V(:);
        K = max(K, 1e-15);
        P = P * sx;

        if abs(sx - 1) < tol * 100, break; end
    end
end

% Newton in log space: u = [ln(K_1),...,ln(K_n), ln(P)]
%
% Residuals:
%   g_i     = ln(K_i) + ln(phi_V_i) - ln(phi_L_i) = 0    (equilibrium)
%   g_{n+1} = sum(z_i / K_i) - 1 = 0                      (summation)

u = [log(max(K, 1e-30)); log(max(P, 1))];
err = Inf;

for loop = 1:maxiter
    g = eval_g(u, comp_vap, temp, pressc, tempc, acentric, BIP, eos_type, components);

    err = max(abs(g));
    if err < tol, break; end

    J = eval_J(u, g, comp_vap, temp, pressc, tempc, acentric, BIP, eos_type, components);

    % Regularize if near-singular
    rcnd = rcond(J);
    if isnan(rcnd) || rcnd < 1e-16
        J = J + 1e-8 * eye(ncomp + 1);
    end

    du = -J \ g;

    % Clamp maximum step in log-space
    max_step = 3.0;
    scale = max(abs(du)) / max_step;
    if scale > 1, du = du / scale; end

    % Backtracking line search
    alpha = 1.0;
    for ls = 1:8
        u_try = u + alpha * du;
        g_try = eval_g(u_try, comp_vap, temp, pressc, tempc, acentric, BIP, eos_type, components);
        if max(abs(g_try)) < err * (1 - 1e-4 * alpha)
            break;
        end
        alpha = alpha * 0.5;
    end

    u = u + alpha * du;
end

K_out = exp(u(1:ncomp));
pressd = exp(u(ncomp + 1));
comp_liq = comp_vap ./ K_out;
comp_liq = max(comp_liq, 0);
comp_liq = comp_liq / sum(comp_liq);

if err > tol
    pressd = NaN;
end

end


function g = eval_g(u, z, T, Pc, Tc, w, BIP, eos, comp)
    ncomp = length(z);
    K = exp(u(1:ncomp));
    P = exp(u(ncomp + 1));

    x = z ./ K;
    x = max(x, 1e-20);
    x = x / sum(x);

    [phi_V, ~] = fugacitycoef_multicomp_vapor(z, P, T, Pc, Tc, w, BIP, eos, comp);
    [phi_L, ~] = fugacitycoef_multicomp_liquid(x, P, T, Pc, Tc, w, BIP, eos, comp);

    g = zeros(ncomp + 1, 1);
    for i = 1:ncomp
        g(i) = u(i) + log(max(phi_V(i), 1e-30)) - log(max(phi_L(i), 1e-30));
    end
    g(ncomp + 1) = sum(z ./ K) - 1;
end


function J = eval_J(u, g0, z, T, Pc, Tc, w, BIP, eos, comp)
    N = length(u);
    J = zeros(N, N);
    for j = 1:N
        h = max(abs(u(j)) * 1e-6, 1e-8);
        u_p = u;
        u_p(j) = u_p(j) + h;
        g1 = eval_g(u_p, z, T, Pc, Tc, w, BIP, eos, comp);
        J(:, j) = (g1 - g0) / h;
    end
end