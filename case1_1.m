clear; clc; close all;

R = 8.3144598;

components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7+'};
comp = [0.0141;0.0014; 0.4930; 0.0978; 0.0467; 0.0062; 0.0132; 0.0022; 0.0296; 0.0182;0.27760]; 
Tc = [126.19 ; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.40;698.985];
Pc = [33.9439e5 ; 73.7646e5; 46.00e5; 48.838e5; 42.455e5; 36.477e5; 37.99e5; 33.8426e5 ; 33.7412e5; 29.6882e5; 19.6872e5];
acentric = [ 0.04 ;0.2250; 0.008; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.296;0.8245];
M_gmol = [ 28.014;44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178;189.393];
ncomp = length(comp);

BIP = zeros(ncomp);

BIP(1,2)= -0.0315; BIP(2,1)= -0.0315;
BIP(1,3)= 0.0278 ; BIP(3,1)= 0.0278 ;
BIP(1,4)= 0.0407 ; BIP(4,1)= 0.0407 ;
BIP(1,5)= 0.0763 ; BIP(5,1)= 0.0763 ;
BIP(1,6)= 0.0944 ; BIP(6,1)= 0.0944 ;
BIP(1,7)= 0.0700 ; BIP(7,1)= 0.0700 ;
BIP(1,8)= 0.0867 ; BIP(8,1)= 0.0867 ;
BIP(1,9)= 0.0878 ; BIP(9,1)= 0.0878 ;
BIP(1,10)= 0.08 ;  BIP(10,1)= 0.08 ;
BIP(1,11)= 0.08 ;  BIP(11,1)= 0.08 ;


for j = 3:ncomp
    BIP(2,j) = 0.12; BIP(j,2) = 0.12;
end

Zc_all = 0.29056 - 0.08775 * acentric;

vt_params.Zc = Zc_all;
vt_params.components = components;

%% Experimental density data - Regueira et al. (2020), Table 3
%  HPHT Volatile Oil, density in g/cm3

T_exp = [298.15; 323.15; 348.15; 358.15; 373.15; 423.15; 432.15; 463.15];
P_exp_MPa = [40; 60; 80; 100; 120; 140];
nT = length(T_exp);
nP = length(P_exp_MPa);

rho_exp = [
    0.6781  0.6945  0.7080  0.7197  0.7301  0.7393
    0.6579  0.6768  0.6918  0.7047  0.7161  0.7261
    0.6374  0.6590  0.6758  0.6899  0.7021  0.7129
    0.6279  0.6505  0.6680  0.6826  0.6950  0.7062
    0.6153  0.6395  0.6581  0.6734  0.6865  0.6980
    0.5755  0.6060  0.6284  0.6463  0.6614  0.6745
    0.5685  0.6002  0.6233  0.6417  0.6572  0.6707
    0.5433  0.5799  0.6056  0.6258  0.6426  0.6570
];

%% Compute model densities

vt_methods = [7, 8, 9];
vt_names = {'SRK (no VT)', 'SRK + Pina-Martinez', 'SRK + Chen-Li'};
nM = length(vt_methods);

rho_model = zeros(nT, nP, nM);

for iT = 1:nT
    for iP = 1:nP
        press = P_exp_MPa(iP) * 1e6;
        Temp  = T_exp(iT);
        for iM = 1:nM
            [rho_kg, ~, ~, ~] = calculate_density(comp, press, Temp, Pc, Tc, acentric, BIP, M_gmol, vt_methods(iM), vt_params);
            rho_model(iT, iP, iM) = rho_kg / 1000;
        end
    end
end

%% Compute AARD for each method

AARD = zeros(nM, 1);
for iM = 1:nM
    rel_err = abs(rho_model(:,:,iM) - rho_exp) ./ rho_exp;
    AARD(iM) = mean(rel_err(:)) * 100;
end

fprintf('%-25s  AARD [%%]\n', 'Method');
fprintf('%s\n', repmat('-', 1, 40));
for iM = 1:nM
    fprintf('%-25s  %.2f\n', vt_names{iM}, AARD(iM));
end

%% Plot: 4 temperature subplots (2x2)

T_plot = [298.15, 323.15, 348.15, 463.15];
T_idx  = [1, 2, 3, 8];

colors = [0.8500 0.3250 0.0980;   % SRK no VT
          0.0000 0.4470 0.7410;   % Pina-Martinez
          0.4660 0.6740 0.1880];  % Chen-Li

line_styles = {'--', '-.', '-'};

figure('Position', [100 100 1200 900]);

for k = 1:4
    subplot(2,2,k);
    hold on; grid on; box on;

    plot(P_exp_MPa, rho_exp(T_idx(k), :), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental');

    for iM = 1:nM
        P_fine = linspace(P_exp_MPa(1), P_exp_MPa(end), 50);
        rho_fine = zeros(1, 50);
        for ip = 1:50
            [rho_kg, ~, ~, ~] = calculate_density(comp, P_fine(ip)*1e6, T_plot(k), Pc, Tc, acentric, BIP, M_gmol, vt_methods(iM), vt_params);
            rho_fine(ip) = rho_kg / 1000;
        end
        plot(P_fine, rho_fine, line_styles{iM}, 'Color', colors(iM,:), 'LineWidth', 2, 'DisplayName', vt_names{iM});
    end

    xlabel('Pressure [MPa]');
    ylabel('\rho [g \cdot cm^{-3}]');
    title(sprintf('T = %.2f K', T_plot(k)));
    if k == 1
        legend('Location', 'southeast', 'FontSize', 7);
    end
    set(gca, 'FontSize', 11);
end

sgtitle('HPHT Volatile Oil - SRK Density Predictions vs Regueira et al. (2020)', 'FontSize', 14, 'FontWeight', 'bold');