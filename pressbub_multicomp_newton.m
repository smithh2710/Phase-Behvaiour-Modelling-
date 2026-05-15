function [pressb, comp_vap] = pressbub_multicomp_newton(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type, components)

if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end
if nargin < 11
    components = {};
end

ncomp = size(comp_liq,1);

pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);
K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);

fun = @(m) objfun(m(1:ncomp, 1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP, eos_type, components);
x = [K; pressb];

for loop = 1:maxiter
    
    f = fun(x);
    J = jacob(f, x, fun);
    dx = -J\f;
    
    x = x + dx;

    err = max(abs(f));
    if err < tol
        break;
    end
    
end

if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp_newton() did not converge: err = %e\n', err);
else
    fprintf('iter = %d, objfun = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

K = x(1:ncomp, 1);
pressb = x(ncomp + 1, 1);

comp_vap = calccompvap(K, comp_liq);

end


function comp_vap = calccompvap(K, comp_liq)

ncomp = size(comp_liq, 1);

comp_vap = zeros(ncomp, 1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

end


function K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components)

ncomp = size(K,1);

comp_vap = calccompvap(K, comp_liq);

[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);

K = zeros(ncomp, 1);
for i = 1:ncomp
    K(i) = fugcoef_liq(i)/fugcoef_vap(i);
end

end


function f = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components)

ncomp = size(K,1);

comp_vap = calccompvap(K, comp_liq);

[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);

f = zeros(ncomp + 1, 1);
for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_vap) - 1;

end


function J = jacob(f0, x, fun)
N = size(x, 1);
J = zeros(N, N);
for i = 1:N
    dx = zeros(N, 1);
    h = max(abs(x(i)) * 1e-7, 1e-10);
    dx(i) = h;
    f1 = fun(x + dx);
    J(:, i) = (f1 - f0) / h;
end
end