function [s, c, mix] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, Zc, components, BIP, verbose)
% CHEN_LI_VOLUME_SHIFT - Volume translation for SRK-EOS using Chen & Li (2020)
% Reference: Chen & Li (2020), Fluid Phase Equilibria 521, 112724
% Alpha function: Pina-Martinez Twu91 with component-specific L, M, N
%
% Uses get_tc_srk_parameters() for L, M, N lookup

R = 8.3144598;

nc = length(comp);
comp = comp(:);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Zc = Zc(:);

if nargin < 8 || isempty(components)
    components = cell(nc, 1);
end

if nargin < 9 || isempty(BIP)
    BIP = zeros(nc);
end

if nargin < 10 || isempty(verbose)
    verbose = false;
end

Zc_SRK = 1/3;
Omega_a = 0.42748;
Omega_b = 0.08664;

table_c123 = {
    'C1',            -0.00195,  0.79540,  2.13497;
    'C2',             0.00270,  0.85431,  2.59463;
    'C3',             0.00492,  0.89221,  2.75570;
    'nC4',            0.00594,  0.93036,  2.60453;
    'iC4',            0.00554,  0.89812,  2.61896;
    'nC5',            0.00831,  0.99529,  3.00384;
    'iC5',            0.00688,  0.89542,  3.07548;
    'nC6',            0.00923,  0.86394,  2.07261;
    'nC7',            0.01032,  1.18249,  2.08864;
    'nC8',            0.01226,  1.19624,  2.35762;
    'nC9',            0.01444,  1.12707,  2.66252;
    'nC10',           0.01574,  1.25347,  3.03108;
    'Cyclohexane',    0.00459,  1.09262,  2.21252;
    'Benzene',        0.00680,  1.00489,  2.55901;
    'Toluene',        0.01024,  1.05282,  2.87184;
    'CO2',            0.00608,  0.92912,  2.65917;
    'N2',            -0.00252,  0.75199,  2.19566;
    'H2S',            0.00144,  0.97009,  2.45887;
    'H2O',            0.02425,  1.30564,  2.17549;
    'NH3',            0.02004,  1.14567,  2.55131;
    'nC11'            0.014409, 1.148682, 3.036902;
    'nC12'            0.017950, 1.064340, 3.480780;
    'nC13'            0.020240, 1.214360, 3.364280;
    'nC14'            0.021119, 1.250016, 2.639114;
    'nc15'            0.023340, 1.194810, 2.886780;
    'nC16'            0.020710, 1.193030, 3.438510;
    'nC17'            0.023390, 1.251400, 3.538670;
    'nC18'            0.023197, 1.119175, 1.697162;
    'nC19'            0.023795, 1.198622, 2.295580;
    'nC20'            0.029300, 1.388140, 3.852120;
    'nC21'            0.027730, 1.290584, 2.586497;
    'nC22'            0.028322, 1.317791, 2.533505;
    'nC23'            0.032306, 1.438407, 3.261193;
    'nC24'            0.031415, 1.314961, 2.258003;
    'nC26'            0.034095, 1.265085, 1.903566;
    'nC27'            0.034095, 1.435975, 3.185357;
    'nC28'            0.032071, 1.352271, 2.531994;   
};


L = zeros(nc, 1);
M = zeros(nc, 1);
N = zeros(nc, 1);
source_LMN = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        [params, found] = get_tc_srk_parameters(components{i});
        if found
            L(i) = params.L;
            M(i) = params.M;
            N(i) = params.N;
            source_LMN{i} = 'Table';
        end
    end
    
    if ~found
        omega_i = acentric(i);
        L(i) = 0.1359 + 0.7535 * omega_i + 0.0611 * omega_i^2;
        M(i) = 0.8787 - 0.2063 * omega_i + 0.1709 * omega_i^2;
        N(i) = 2.0;
        source_LMN{i} = 'Generalized';
    end
end

%% Get c1, c2, c3 for each component
c1 = zeros(nc, 1);
c2 = zeros(nc, 1);
c3 = zeros(nc, 1);
source_c = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        comp_name = strtrim(components{i});
        
        for j = 1:size(table_c123, 1)
            if strcmpi(comp_name, table_c123{j, 1})
                c1(i) = table_c123{j, 2};
                c2(i) = table_c123{j, 3};
                c3(i) = table_c123{j, 4};
                source_c{i} = 'Table';
                found = true;
                break;
            end
        end
    end
    
    if ~found
        omega_i = acentric(i);
        Zc_i = Zc(i);
        
        c1(i) = -3.90812e-4 + 0.03274 * omega_i;
        
        if Zc_i < 0.2215
            Zc_i = 0.2215;
        elseif Zc_i > 0.2864
            Zc_i = 0.2864;
        end
        
        c2(i) = 3.06048 - 7.64314 * Zc_i;
        c3(i) = 8.34576 - 21.07619 * Zc_i;
        
        source_c{i} = 'Generalized';
    end
end

%% Calculate EOS parameters with component-specific Twu91 alpha
Tr = temp ./ Tc;
alpha = Tr.^(N.*(M-1)) .* exp(L .* (1 - Tr.^(N.*M)));

a_i = Omega_a * R^2 * Tc.^2 ./ Pc .* alpha;
b_i = Omega_b * R * Tc ./ Pc;

a_mix = 0;
for i = 1:nc
    for j = 1:nc
        a_mix = a_mix + comp(i) * comp(j) * sqrt(a_i(i) * a_i(j)) * (1 - BIP(i,j));
    end
end
b_mix = sum(comp .* b_i);

%% Solve SRK cubic for Z
A = a_mix * press / (R^2 * temp^2);
B = b_mix * press / (R * temp);

coeffs = [1, -1, (A - B - B^2), -A*B];
Z_roots = roots(coeffs);
Z_roots = Z_roots(imag(Z_roots) == 0 & real(Z_roots) > 0);

if isempty(Z_roots)
    Z = 0.3;
else
    Z = min(real(Z_roots));
end

V_EOS = Z * R * temp / press;

%% Fugacity coefficients from SRK with Twu91 alpha (self-consistent)
fugcoef = zeros(nc, 1);
for i = 1:nc
    sum_xaij = 0;
    for j = 1:nc
        sum_xaij = sum_xaij + comp(j) * sqrt(a_i(i) * a_i(j)) * (1 - BIP(i,j));
    end
    ln_phi = (b_i(i)/b_mix)*(Z-1) - log(Z-B) ...
        - (A/B)*(2*sum_xaij/a_mix - b_i(i)/b_mix)*log(1 + B/Z);
    fugcoef(i) = exp(ln_phi);
end

%% For pure component or mixture
if nc == 1
    Tc_use = Tc(1);
    Pc_use = Pc(1);
    Zc_exp = Zc(1);
    seta = 1;
    Vc_mix = Zc(1) * R * Tc(1) / Pc(1);
else
    Vc_est = Zc .* R .* Tc ./ Pc;
    Vc_23 = Vc_est.^(2/3);
    sum_xVc23 = sum(comp .* Vc_23);
    if sum_xVc23 > 1e-15
        seta = (comp .* Vc_23) / sum_xVc23;
    else
        seta = comp;
    end
    Tc_use = sum(seta .* Tc);
    Vc_mix = sum(seta .* Vc_est);
    omega_mix = sum(comp .* acentric);
    Pc_use = (0.2905 - 0.085 * omega_mix) * R * Tc_use / Vc_mix;
    Zc_exp = 0.2905 - 0.085 * omega_mix;
end

%% Distance function (Eq. 5)
term1 = (V_EOS^2 * temp) / (Tc_use * (V_EOS - b_mix)^2);
term2 = (a_mix * (2*V_EOS + b_mix)) / (R * Tc_use * (V_EOS + b_mix)^2);
d_mix = term1 - term2;

%% Volume translation parameters
c1_mix = sum(comp .* c1);
c2_mix = sum(comp .* c2);
c3_mix = sum(comp .* c3);

%% delta_c (Eq. 7)
delta_c = (R * Tc_use / Pc_use) * (Zc_SRK - Zc_exp);

%% Total volume translation (Eq. 9)
c_far = c1_mix * R * Tc_use / Pc_use;
c_near = delta_c / (c2_mix + c3_mix * d_mix);
c_total = c_far + c_near;

V_VTEOS = V_EOS - c_total;

%% Output
c = c1 .* R .* Tc ./ Pc;
s = c ./ b_i;

mix.c_mix = c_total;
mix.c_mix_cm3 = c_total * 1e6;
mix.delta_c = delta_c;
mix.d = d_mix;
mix.c1 = c1;
mix.c2 = c2;
mix.c3 = c3;
mix.c1m = c1_mix;
mix.c2m = c2_mix;
mix.c3m = c3_mix;
mix.L = L;
mix.M = M;
mix.N = N;
mix.alpha = alpha;
mix.Tc_mix = Tc_use;
mix.Pc_mix = Pc_use;
mix.Vc_mix = Vc_mix;
mix.V_EOS = V_EOS;
mix.V_VTEOS = V_VTEOS;
mix.Z = Z;
mix.source_LMN = source_LMN;
mix.source_c = source_c;
mix.seta = seta;
mix.a = a_mix;
mix.b = b_mix;
mix.fugcoef = fugcoef;

if verbose
    fprintf('Chen-Li Volume Shift Results:\n');
    fprintf('  V_EOS = %.6e m3/mol\n', V_EOS);
    fprintf('  c_total = %.6e m3/mol (%.4f cm3/mol)\n', c_total, c_total*1e6);
    fprintf('  V_VTEOS = %.6e m3/mol\n', V_VTEOS);
    for i = 1:nc
        fprintf('  Comp %d: L=%.4f, M=%.4f, N=%.4f (%s)\n', i, L(i), M(i), N(i), source_LMN{i});
        fprintf('          c1=%.5f, c2=%.5f, c3=%.5f (%s)\n', c1(i), c2(i), c3(i), source_c{i});
    end
end

end