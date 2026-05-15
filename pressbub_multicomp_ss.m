function [pressb, comp_vap] = pressbub_multicomp_ss(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type, components)

if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end
if nargin < 11
    components = {};
end

ncomp = size(comp_liq,1);

pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);

for loop = 1:maxiter
    
    [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);
    
    K = Knew;
    pressb = pressbnew;
    
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp_ss() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d\n', loop);
end

comp_vap = calccompvap(K, comp_liq);

end


function comp_vap = calccompvap(K, comp_liq)

ncomp = size(comp_liq, 1);

comp_vap = zeros(ncomp, 1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

end


function [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components)

ncomp = size(K,1);

comp_vap = calccompvap(K, comp_liq);

[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type, components);

fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);
pressbnew = 0;
for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressb;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressb;
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    pressbnew = pressbnew + fug_liq(i)/fugcoef_vap(i);
    
end

f = comp_liq'*Knew - 1;

end