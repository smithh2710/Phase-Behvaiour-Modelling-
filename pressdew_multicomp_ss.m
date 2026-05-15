function [pressd, comp_liq, K_out] = pressdew_multicomp_ss(comp_vap, pressd_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type, components, K_init)

if nargin < 10 || isempty(eos_type), eos_type = 'PR'; end
if nargin < 11, components = {}; end

ncomp = size(comp_vap, 1);
comp_vap = comp_vap(:);

have_K_init = (nargin >= 12) && ~isempty(K_init) && all(K_init > 0) && all(isfinite(K_init));

if have_K_init
    K = K_init(:);
else
    K = wilsoneq(pressd_ini, temp, pressc, tempc, acentric);
end

pressd = pressd_ini;
eps_val = Inf;

for loop = 1:maxiter

    [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type, components);

    K = Knew;
    pressd = pressdnew;

    eps_val = abs(f);
    if eps_val < tol
        break;
    end

end

K_out = K;

if loop >= maxiter
    pressd = NaN;
end

comp_liq = comp_vap ./ K;
comp_liq = max(comp_liq, 0);
comp_liq = comp_liq / sum(comp_liq);

end


function comp_liq = calccompliq(K, comp_vap)
    ncomp = size(comp_vap, 1);
    comp_liq = zeros(ncomp, 1);
    for i = 1:ncomp
        comp_liq(i) = comp_vap(i) / K(i);
    end
end


function [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type, components)

    ncomp = size(K, 1);

    comp_liq = calccompliq(K, comp_vap);
    comp_liq = max(comp_liq, 1e-20);
    comp_liq = comp_liq / sum(comp_liq);

    [fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type, components);
    [fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP, eos_type, components);

    Knew = zeros(ncomp, 1);
    f = 0;
    pressdnew = 0;

    for i = 1:ncomp
        fug_vap_i = comp_vap(i) * fugcoef_vap(i) * pressd;
        Knew(i) = fugcoef_liq(i) / fugcoef_vap(i);
        Knew(i) = max(Knew(i), 1e-15);
        pressdnew = pressdnew + fug_vap_i / fugcoef_liq(i);
        f = f + comp_vap(i) / Knew(i);
    end

    f = f - 1;

end