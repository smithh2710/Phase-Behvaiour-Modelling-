function tau = esposito_tau_correlation(comp, T, M_gmol, Tc, Pc, params)

alpha0 = params(1);
alpha1 = params(2);
beta0  = params(3);
beta1  = params(4);

comp = comp(:);
comp = comp / sum(comp);
M_gmol = M_gmol(:);
Tc = Tc(:);
Pc = Pc(:);

Tc_mix = sum(comp .* Tc);
Tr_mix = T / Tc_mix;

alpha = alpha0 + alpha1 * Tr_mix;
beta  = beta0  + beta1  / Tr_mix;

ncomp = length(M_gmol);
tau = zeros(ncomp, 1);
for i = 1:ncomp
    tau(i) = M_gmol(i)^beta * exp(alpha);
end

tau = max(tau, 0.5);
tau = min(tau, 10.0);

end