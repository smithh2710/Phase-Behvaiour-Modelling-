components = {'CO2','N2','C1','C2','C3','iC4','nC4','iC5','nC5','nC6','nC11','nC33'};
comp = [0.12 , 0.13 , 54.98 , 5.44 , 4.88 ,1.88 , 2.48 , 1.45 , 1.30 , 1.91 , 21.53 , 3.90 ]' ; 
comp = comp./sum(comp) ; 
pressc = [73.76 , 33.94 , 46 , 48.84 , 42.46 , 36.48 , 38 , 33.84 , 33.74 , 29.69 , 23.43 , 9.90 ]'*10^5 ; 
tempc = [304.2 , 126.2 , 190.60 , 305.4 , 369.8 , 408.1 , 425.2 , 460.4 , 469.6 , 507.4 , 642.51 , 929.34]'; 
accentric = [0.2250 , 0.04 , 0.008 , 0.098 , 0.1520 , 0.1760 , 0.1930 , 0.227 , 0.2510 , 0.2960 , 0.4882 , 1.0854 ]' ; 
M_gmol = [44.01 , 28.01 , 16.04 , 30.07 , 44.1, 58.12 , 58.12 ,72.15 , 72.15 , 86.18 , 151.16 , 459.34]' ; 

BIP = zeros(n);
BIP(3,4:end) = [0.0256, 0.0270, 0.0284, 0.0284, 0.0298, 0.0298, 0.0312, 0.0486, 0.083];
BIP(4:end,3) = BIP(3,4:end);

h_ref = 100;
press_ref = 293.02e5;
temp_ref = 361;
dTdh = 0;
eos_type = 'PR';
tol = 1e-10; 
maxiter = 500;

Cp_coeffs = zeros(n, 4);
H_ig_ref = zeros(n, 1);

depth_range = [0, 120];

vt_methods = [0];
vt_labels  = {'PR (no VT)', 'Baled', 'Abudour'};
n_vt = length(vt_methods);

Zc = 0.29056 - 0.08775 .* acentric;
vt_params = struct();
vt_params.components = components;
vt_params.Zc = Zc;


[GOC_depth, GOC_pressure, ~, GOC_temp] = detect_sgoc(1, comp, press_ref, temp_ref,h_ref, dTdh, depth_range, pressc, tempc, accentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, 0, vt_params)