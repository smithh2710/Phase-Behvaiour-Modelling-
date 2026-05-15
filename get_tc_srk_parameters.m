function [params, found] = get_tc_srk_parameters(comp_name)
% GET_TC_SRK_PARAMETERS - Lookup table for tc-SRK (Pina-Martinez 2018)
%
% Returns L, M, N (Twu91 alpha) and c (volume translation) for known components
%
% USAGE:
%   [params, found] = get_tc_srk_parameters('H2O')
%   [params, found] = get_tc_srk_parameters('C1')
%
% OUTPUT:
%   params.L    - Twu91 L parameter
%   params.M    - Twu91 M parameter
%   params.N    - Twu91 N parameter
%   params.c    - Volume translation [cm3/mol]
%   params.name - Matched component name
%   found       - true if component found in table
%
% REFERENCES:
%   Pina-Martinez et al. (2018) J. Chem. Eng. Data 63, 3980-3988
%   Table S2: tc-RK (SRK) parameters

% tc-SRK parameters from Pina-Martinez (2018) Supporting Information Table S2
% Format: {Name, L, M, N, c [cm3/mol]}
table = {
    'N2',             0.1901,  0.8900,  2.0107,   1.3475;
    'CO2',            0.2806,  0.8684,  2.2782,   4.1585;
    'H2S',            0.1748,  0.8686,  2.2761,   3.0181;
    'H2O',            0.4171,  0.8758,  2.1818,   8.9670;
    'C1',             0.2170,  0.9082,  1.8172,   2.0509;
    'C2',             0.2968,  0.8812,  1.7252,   4.6079;
    'C3',             0.5427,  0.8811,  1.1904,   7.7000;
    'iC4',            0.6853,  0.8842,  1.0305,  10.8558;
    'nC4',            0.3515,  0.8609,  1.8323,  11.0274;
    'neoC5',          0.1821,  0.8476,  2.5887,  11.6811;
    'iC5',            0.2507,  0.8548,  2.3951,  13.8762;
    'nC5',            0.2950,  0.8513,  2.2388,  16.2371;
    '22DMC4',         0.2316,  0.8523,  2.5126,  13.9389;
    '23DMC4',         0.2391,  0.8510,  2.5337,  16.1281;
    '2MC5',           0.3279,  0.8404,  2.1294,  17.7302;
    '3MC5',           0.3141,  0.8389,  2.1441,  19.0017;
    'nC6',            0.2982,  0.8491,  2.4179,  22.0561;
    'nC7',            0.3529,  0.8368,  2.2644,  27.8526;
    'nC8',            0.3597,  0.8321,  2.3840,  34.8744;
    'nC9',            0.4046,  0.8255,  2.3023,  41.4328;
    'nC10',           0.3786,  0.8255,  2.5740,  48.7183;
    'nC11',           0.5303,  0.8488,  2.2228,  58.4006;
    'nC12',           0.4299,  0.8171,  2.5617,  65.9614;
    'nC13',           0.4659,  0.8143,  2.5169,  78.4645;
    'nC14',           0.5090,  0.8079,  2.3926,  92.7429;
    'nC15',           0.5292,  0.8044,  2.4124, 102.3182;
    'nC16',           0.5435,  0.8042,  2.4387, 113.3416;
    'nC17',           0.5331,  0.8075,  2.6166, 117.2450;
    'nC18',           0.5664,  0.8054,  2.5832, 128.1272;
    'nC19',           0.5824,  0.8026,  2.6145, 137.9881;
    'nC20',           0.5596,  0.8140,  2.8872, 140.6330;
    'nC21',           0.5191,  0.8173,  3.1709, 153.5489;
    'nC23',           0.5437,  0.8142,  3.2430, 173.8344;
    'nC24',           0.5559,  0.8124,  3.2859, 183.0294;
    'nC26',           0.6098,  0.8135,  3.2596, 201.1113;
    'nC27',           0.6766,  0.8193,  3.1268, 205.7583;
    'nC28',           0.6633,  0.8149,  3.2252, 217.2250;
    'Cyclopentane',   0.3229,  0.8496,  1.8673,  12.6704;
    'Cyclohexane',    0.3133,  0.8336,  1.8785,  13.3342;
    'Benzene',        0.1919,  0.8469,  2.5997,  13.4541;
    'Toluene',        0.3312,  0.8448,  2.0808,  19.6570;
};

params = struct('L', [], 'M', [], 'N', [], 'c', [], 'name', '');
found = false;

if isempty(comp_name)
    return;
end

comp_name = strtrim(comp_name);

for i = 1:size(table, 1)
    if strcmpi(comp_name, table{i, 1})
        params.L = table{i, 2};
        params.M = table{i, 3};
        params.N = table{i, 4};
        params.c = table{i, 5};
        params.name = table{i, 1};
        found = true;
        return;
    end
end

end