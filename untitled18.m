components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','nC6','nC7','nC8','nC9','nC10','nC14','nC17','nC20','nC23','nC28','nC34','nC39-48','C49-C80'};
      
comp_ref = [0.395; 2.060; 53.8715362; 7.589; 5.575; 1.009; 2.514; 0.900; 1.396; 1.557;2.630; 2.823; 1.897; 4.406; 2.479; 1.941; 1.520; 1.526; 1.325; 1.145;0.829; 0.610] ;
comp_ref = comp_ref./(sum(comp_ref)); 
M_gmol = [28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178;96.000; 107.000; 121.000; 152.857; 205.131; 249.627; 289.518; 336.857; 399.72320; 481.458;595.685; 811.198];


Tc = [-146.950; 31.050; -82.550; 32.250; 96.650; 134.950; 152.050; 187.250; 196.450; 234.250;262.184; 282.716; 304.023; 346.593; 401.967; 442.808; 476.719; 514.042; 560.087; 616.024;689.007; 821.700] + 273.15;
Pc = [33.958; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; 29.69;31.95; 29.76; 26.67; 22.04; 18.06; 16.32; 15.43; 14.71; 14.12; 13.69;13.43; 13.40] * 1e5;
acentric = [0.0377; 0.2250; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.2960;0.4679; 0.4999; 0.5399; 0.6321; 0.7673; 0.8746; 0.9642; 1.0608; 1.1710; 1.2799;1.3595; 1.2257];


c_custom = [0.92; 3.03; 0.63; 2.63; 5.06; 7.29; 7.86; 10.93; 12.18; 17.98;9.65; 13.90; 20.99; 34.32; 45.88; 45.87; 39.23; 26.41; 3.39; -33.93;-94.44; -219.18];
 Vc = [89.80; 94.00; 99.00; 148.00; 203.00; 263.00; 255.00; 306.00; 304.00; 370.00;432; 487.84; 529.18; 603; 830; 1000; 1140; 1305; 1580; 1910; 2250; 4022.43]*10^(-6);
 Zc = Pc .* (Vc ) ./ (8.3145 * Tc);

%  Zc_exp = [0.2894; 0.2744; 0.2862; 0.2793; 0.2763; 0.2824; 0.2739; 0.2701; 0.2685; 0.2635;0.2611; 0.2559; 0.2520];
% Zc_pseudo = 0.29056 - 0.08775.*acentric(14:22);
% for i = 2:length(Zc_pseudo)
%     if Zc_pseudo(i) >= Zc_pseudo(i-1)
%         Zc_pseudo(i) = Zc_pseudo(i-1) - 0.008;
%     end
% end
% Zc = [Zc_exp ; Zc_pseudo];
R = 8.314456 ;
% Vc = Zc .* R .* Tc ./ Pc;
vt_params.c_custom = c_custom;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

% Tc = [126.200; 304.210; 190.564; 305.320; 369.830; 407.800; 425.120; 460.400; 469.700; 507.600; ...
%       262.184+273.15; 282.716+273.15; 304.023+273.15; 346.593+273.15; 401.967+273.15; 442.808+273.15; ...
%       476.719+273.15; 514.042+273.15; 560.087+273.15; 616.024+273.15; 689.007+273.15; 821.700+273.15];
% 
% Pc = [34.000; 73.830; 45.990; 48.720; 42.480; 36.400; 37.960; 33.800; 33.700; 30.250; ...
%       31.95; 29.76; 26.67; 22.04; 18.06; 16.32; 15.43; 14.71; 14.12; 13.69; 13.43; 13.40] * 1e5;
% 
% acentric = [0.0377; 0.2236; 0.0115; 0.0995; 0.1523; 0.1835; 0.2002; 0.2279; 0.2515; 0.3013; ...
%             0.4679; 0.4999; 0.5399; 0.6321; 0.7673; 0.8746; 0.9642; 1.0608; 1.1710; 1.2799; 1.3595; 1.2257];



n = 22;
h_ref    = 3682.8;   
P_ref    = 380.20e5;     
T_ref    = 138.70 + 273.15; 
dTdh     = 0.0260; 

BIP = zeros(n, n);
BIP(2,1)  = -0.0315; 
BIP(3,1)  =  0.0278;
BIP(4,1)  =  0.0407;
BIP(5,1)  =  0.0763;
BIP(6,1)  =  0.0944;
BIP(7,1)  =  0.0700;  
BIP(8,1)  =  0.0867;  
BIP(9,1)  =  0.0878;  
BIP(10,1) =  0.0800; 

for i = 11:22
    BIP(i,1) = 0.0800; 
end

for i = 3:10
    BIP(i,2) = 0.1200;  
end
for i = 11:22
    BIP(i,2) = 0.1000;
end

BIP = BIP + BIP';


Cp_coeffs = [
    31.15,   -0.014,    2.68e-5,  -1.17e-8;   % N2
    19.79,    0.073,   -5.60e-5,   1.72e-8;   % CO2
    19.25,    0.052,    1.20e-5,  -1.13e-8;   % C1
     5.41,    0.178,   -6.94e-5,   8.71e-9;   % C2
    -4.22,    0.306,   -1.59e-4,   3.21e-8;   % C3
    -1.39,    0.385,   -1.85e-4,   2.90e-8;   % iC4
     9.49,    0.331,   -1.11e-4,  -2.82e-9;   % nC4
     9.52,    0.507,   -2.73e-4,   5.72e-8;   % iC5
    -3.63,    0.487,   -2.58e-4,   5.30e-8;   % nC5
    -4.41,    0.582,   -3.12e-4,   6.49e-8;   % C6
     9.58,    0.576,   -2.05e-4,   0.000;     % C7
    -4.91,    0.612,   -2.33e-4,   0.000;     % C8
    -1.87,    0.680,   -2.69e-4,   0.000;     % C9
     2.31,    0.859,   -3.48e-4,   0.000;     % C10-C13
     5.06,    1.151,   -4.67e-4,   0.000;     % C14-C16
     7.19,    1.402,   -5.66e-4,   0.000;     % C17-C19
     8.41,    1.629,   -6.57e-4,   0.000;     % C20-C22
    10.06,    1.899,   -7.65e-4,   0.000;     % C23-C26
    12.3129,    2.258,   -9.08e-4,   0.000;     % C27-C31
    14.70,    2.721,   -1.09e-3,   0.000;     % C32-C38
    18.27,    3.369,   -1.36e-3,   0.000;     % C39-C48
    26.84,    4.652,   -1.87e-3,   0.000];    % C49-C80


H_ig_ref = [8330.8; 19459.1; 2.6; 9761.1; 19519.6; 29278.1; 29278.1; 39036.6; 39036.6; 48795.1;55628.2; 63280.8; 73020.5; 95183.2; 131550.0; 162505.2; 190257.4; 223190.8; 266926.4617737;323789.1; 403255.9; 553186.7];

% Zc_exp = [0.2894; 0.2744; 0.2862; 0.2793; 0.2763; 0.2824; 0.2739; 0.2701; 0.2685; 0.2635;0.2611; 0.2559; 0.2520];

% Zc = [0.2894, 0.2744, 0.2862, 0.2793, 0.2763, 0.2824, 0.2739, 0.2716, 0.2685, 0.2635, 0.2611, 0.2559, 0.2520, 0.2494, 0.2368, 0.2285, 0.2215, 0.2150, 0.2050, 0.1940, 0.1850, 0.1700];

% Zc_pseudo = Pc(14:22) .* (Vc(14:22) * 1e-6) ./ (8.3145 * Tc(14:22));
% Zc = [Zc_exp; Zc_pseudo];

vt_params.c_custom = c_custom;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;

BIP = zeros(n, n);

BIP(2,1)  = -0.0315; 
BIP(3,1)  =  0.0278;
BIP(4,1)  =  0.0407;
BIP(5,1)  =  0.0763;
BIP(6,1)  =  0.0944;
BIP(7,1)  =  0.0700;  
BIP(8,1)  =  0.0867;  
BIP(9,1)  =  0.0878;  
BIP(10,1) =  0.0800; 

for i = 11:22
    BIP(i,1) = 0.0800; 
end

for i = 3:10
    BIP(i,2) = 0.1200;  
end
for i = 11:22
    BIP(i,2) = 0.1000;
end

BIP = BIP + BIP';

h = [3638.2;  3644.3; 3651.1;  3661.6;  3676.0;  3682.8]; 
vt_method = 7 ;

 % [comp_h, press_h, temp_h,~, ~] = main_hasse(h(3) , h(6), comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)

 [GOC_depth, press_h, comp_h,temp_h] = detect_sgoc(3, comp_ref, P_ref, T_ref, h_ref, dTdh, [3600 3660], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)
% [rho, ~, ~, ~] = calculate_density(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, M_gmol, vt_method , vt_params)

% [s, c, mix] = chen_li_volume_shift(comp_ref, P_ref, T_ref, Pc, Tc, acentric, Zc, components)


 % [press_bub, ~] = pressbub_multicomp_ss(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, 1e-6, 500,'SRK', components)

% comp_gas_GOC = comp_gas_GOC *100 

 % [stability, x_trial1, x_trial2, TPD_min, converged] = stability_analysis_ssi(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, 1e-6,500,'SRK' ,components)


%% GOC detection and composition for 3 SRK VT methods
R = 8.3144598;


% c_peneloux = [-4.23; -1.64; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12; 1.39; ...
%                6.70; 10.80; 14.52; 21.55; 26.56; 25.75; 22.30; 15.13; 1.78; -19.81; -55.16; -124.90];

% vt_methods = [7, 12, 9];
vt_methods = [7];
vt_names   = {'SRK (NO VT)', 'SRK + Custom VT', 'SRK + Chen-Li'};

for k = 1:1
    vm = vt_methods(k);

    vp = vt_params;
    if vm == 7
        vp.c_custom = c_custom;
    end


    % [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = detect_sgoc(1, comp_ref, P_ref, T_ref, h_ref, dTdh, [3600 3650], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vm, vp);



    % [~, comp_gas_GOC] = pressbub_multicomp_newton(GOC_comp, GOC_pressure, GOC_temp, Pc, Tc, acentric, BIP, 1e-10, 500)
    % 
    % [~, comp_gas__GOC] = pressbub_multicomp_ss(GOC_comp, GOC_pressure, GOC_temp, Pc, Tc, acentric, BIP, 1e-10, 500)
    % 
    % [stability, x_trial1, x_trial2, TPD_min, converged] = stability_analysis_ssi(GOC_comp, GOC_pressure, GOC_temp, Pc, Tc, acentric, BIP, 1e-4,500,'SRK' ,components)
    % 
    % [stability1, x__trial1, x__trial2, TPD__min] = stability_analysis_qnss(GOC_comp, GOC_pressure, GOC_temp, Pc, Tc, acentric , BIP, 1e-4,500,'SRK' )

    % [K, comp_vap, comp_liq, phasefrac] = vaporliquideq(GOC_pressure, GOC_temp, GOC_comp, Pc, Tc, acentric, BIP, 1e-4, 500,'SRK', components)

end