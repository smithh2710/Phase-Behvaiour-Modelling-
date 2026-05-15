
% charcaterized for SRK EOS , field data from Creek and Scharder 1985 
% Underdsaturated GOC 

components = {'N2' , 'CO2', 'C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7','C8','C9','F1','F2','F3' }; 
comp_ref = [0.131 ; 0.192 ; 59.475 ; 12.010 ; 6.859 ; 1.838 ; 2.485 ; 1.071 ; 1.0 ; 1.404 ; 1.884 ; 1.622 ; 1.396 ;4.553 ; 2.420 ; 1.660 ] ; 
comp_ref = comp_ref./sum (comp_ref) ; 
M_gmol = [28.014 ; 44.010 ; 16.043 ; 30.070 ; 44.097 ; 58.124 ; 58.124 ; 72.151 ; 72.151 ; 86.178 ; 96 ; 107 ; 121 ; 157.252 ; 236.356 ; 376.503 ]; 
Pc = [33.94 ; 73.76 ; 46 ; 48.84 ; 42.46 ; 36.48 ; 38 ; 33.84 ; 33.74 ; 29.69 ; 29.84 ; 27.13 ; 24.27 ;19.94 ; 15.80 ; 13.91]*10^5 ; % Pa
Tc = [126.2 ; 304.2 ; 190.6 ; 305.4 ; 369.8 ;408.1 ; 425.2; 460.4 ; 469.6 ; 507.4 ; 530.873 ; 549.838 ; 571.037 ; 620.361 ; 701.503 ; 826.633 ] ; 
acentric = [0.04 ; 0.225 ; 0.008 ; 0.098 ; 0.152 ; 0.176 ; 0.193 ; 0.2270 ; 0.2510 ; 0.2960 ; 0.4677 ; 0.4996 ; 0.5396 ; 0.6459 ; 0.8465 ; 1.1333 ]; 

ncomp = length(comp_ref); 
BIP = zeros(ncomp,ncomp);
BIP (1,2:ncomp) = [-0.0315 ; 0.0278 ; 0.0407 ; 0.0763 ; 0.0944 ; 0.0700; 0.0867 ; 0.0878 ; 0.08 ; 0.080 ; 0.080; 0.080 ; 0.080 ; 0.080 ; 0.080];
BIP (2,3:ncomp) = [0.12;0.12;0.12;0.12;0.12;0.12;0.12;0.12;0.1;0.1;0.1;0.1;0.1;0.1;];
BIP = BIP + BIP' ; 
press_ref = 317e5 ; % Pa 
temp_ref = 88 + 273.15 ; % K 
h_ref = 1388 ; 
depth_range = [1200 ; 1400  ] ; 
dTdh = 0 ; 
Zc = 0.29056 - 0.08775 .* acentric; 
Cp_coeffs = zeros(ncomp,ncomp) ; 
H_ig_ref = zeros(ncomp);
vt_method = 7 ; 
vt_params.components = components; 
vt_params.Zc = Zc;


[GOC_depth, GOC_pressure, GOC_temp, GOC_comp, GOC_type, info] = detect_GOC(comp_ref, press_ref, temp_ref, h_ref, depth_range, ...
    Pc, Tc, acentric, BIP, M_gmol, ...
    dTdh, Cp_coeffs, H_ig_ref, vt_method, vt_params)