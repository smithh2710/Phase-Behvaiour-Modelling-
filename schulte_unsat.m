
components = {'C1','C2','C3','nC5','nC10'};
comp_ref = [80.97 ; 5.66 ; 3.06 ; 4.57 ; 3.30;2.44 ];
comp_ref = comp_ref ./ sum(comp_ref); 
Pc = [46 ; 48.84 ; 42.46 ;33.74 ; 27.36 ; 21.08]*1e5; 
Tc = [190.6 ; 305.4 ; 369.8 ; 469.6 ; 540.2 ; 617.6 ]; 
acentric = [0.008 ; 0.098 ; 0.152 ; 0.251 ; 0.351 ; 0.49]; 
M_gmol = [16.04 ; 30.07 ; 44.10 ; 72.15 ; 100.21 ; 142.29];

n = length(comp_ref); 

BIP = zeros(n,n);
BIP(1, 2:6) = [0.0256, 0.027, 0.0298, 0.0326, 0.0368];
BIP = BIP + BIP';
Zc = 0.29056 - 0.08775 .* acentric;

h_ref = 1000;
press_ref = 2275e5;
temp_ref = 275;
dTdh = 0;
eos_type = 'PR';
tol = 1e-10;
maxiter = 500;

Cp_coeffs = zeros(n, 4);
H_ig_ref = zeros(n, 1);

depth_range = [900, 1200];
    vt_params = struct();
    vt_params.components = components;
    vt_params.Zc = Zc;
    vt_params.Vc = Zc .* R .* Tc ./ Pc * 1e6;

 [GOC_depth, GOC_pressure, GOC_temp, GOC_comp, GOC_type, info] = detect_GOC(...
    comp_ref, press_ref, temp_ref, h_ref, depth_range, ...
    Pc, Tc, acentric, BIP, M_gmol, ...
    dTdh, Cp_coeffs, H_ig_ref, 0, vt_params)