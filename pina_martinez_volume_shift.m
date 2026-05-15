function [s, c, mix] = pina_martinez_volume_shift(T, Pc, Tc, omega, components, comp)

R = 8.314472;
Omega_b = 0.08664;

nc = length(components);
Pc = Pc(:);
Tc = Tc(:);
omega = omega(:);

if nargin < 6 || isempty(comp)
    comp = ones(nc, 1) / nc;
else
    comp = comp(:);
end

c_cm3 = zeros(nc, 1);
source = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        [params, found] = get_tc_srk_parameters(components{i});
        if found
            c_cm3(i) = params.c;
            source{i} = 'Table';
        end
    end
    
    if ~found
        Zra = 0.29056 - 0.08775 * omega(i);
        c_cm3(i) = 0.40768 * (0.29441 - Zra) * R * Tc(i) / Pc(i) * 1e6;
        source{i} = 'Correlation';
    end
end

c = c_cm3 * 1e-6;

b_i = Omega_b * R * Tc ./ Pc;
s = c ./ b_i;

c_mix = sum(comp .* c);

mix = struct();
mix.c_i = c;
mix.c_i_cm3 = c_cm3;
mix.s = s;
mix.source = source;
mix.c_mix = c_mix;
mix.c_mix_cm3 = c_mix * 1e6;

end