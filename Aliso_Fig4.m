clc; clear; close all;

R = 8.3144598;


h_top    = 2290.3;
h_bot    = 2440;
h_ref    = 2351.3;
temp_ref = 352.32;
dTdh     = 0;
press_ref = 280.00e5;   


comp_ref = [10; 1.917; 1.098; 63.477; 7.821; 5.049; 7.443; 2.781; 0.396; 0.018] ;
comp_ref = comp_ref / sum(comp_ref);   

M_gmol  = [2.016;  28.014; 44.010; 16.043; 30.070; 44.097; 58.199; 107.800; 198.500; 335.100];
Tc      = [33.20;  126.20; 304.18; 190.60; 305.32; 369.83; 460.00; 578.000; 718.000; 868.000];
Pc      = [12.97;  33.94;  73.76;  46.00;  48.84;  42.46;  33.50;  24.800;  16.200;  11.800] * 1e5;
acentric = [-0.220; 0.040;  0.225;  0.008;  0.098;  0.152;  0.1933;  0.412;   0.723;   1.050];

nc = 10;
BIP = zeros(nc);


BIP(1,2) =  0.0069; BIP(2,1) =  0.0069;   % H2-N2
BIP(1,3) = -0.3426; BIP(3,1) = -0.3426;   % H2-CO2
BIP(1,4) = -0.0220; BIP(4,1) = -0.0220;   % H2-C1
BIP(1,5) = -0.1667; BIP(5,1) = -0.1667;   % H2-C2
BIP(1,6) = -0.2359; BIP(6,1) = -0.2359;   % H2-C3
BIP(1,7) = -0.5071; BIP(7,1) = -0.5071;   % H2-iC4nC6
BIP(1,8) =  0.2500; BIP(8,1) =  0.2500;   % H2-F1
BIP(1,9) =  0.2500; BIP(9,1) =  0.2500;   % H2-F2
BIP(1,10)=  0.2500; BIP(10,1)=  0.2500;   % H2-F3
% N2 BIPs (index 2)
BIP(2,3) = -0.0315; BIP(3,2) = -0.0315;   % N2-CO2
BIP(2,4) =  0.0278; BIP(4,2) =  0.0278;   % N2-C1
BIP(2,5) =  0.0407; BIP(5,2) =  0.0407;   % N2-C2
BIP(2,6) =  0.0763; BIP(6,2) =  0.0763;   % N2-C3
BIP(2,7) =  0.0702; BIP(7,2) =  0.0702;   % N2-iC4nC6
BIP(2,8) =  0.0800; BIP(8,2) =  0.0800;   % N2-F1
BIP(2,9) =  0.0800; BIP(9,2) =  0.0800;   % N2-F2
BIP(2,10)=  0.0800; BIP(10,2)=  0.0800;   % N2-F3
% CO2 BIPs (index 3)
BIP(3,4) =  0.1200; BIP(4,3) =  0.1200;   % CO2-C1
BIP(3,5) =  0.1200; BIP(5,3) =  0.1200;   % CO2-C2
BIP(3,6) =  0.1200; BIP(6,3) =  0.1200;   % CO2-C3
BIP(3,7) =  0.1200; BIP(7,3) =  0.1200;   % CO2-iC4nC6
BIP(3,8) =  0.1000; BIP(8,3) =  0.1000;   % CO2-F1
BIP(3,9) =  0.1000; BIP(9,3) =  0.1000;   % CO2-F2
BIP(3,10)=  0.1000; BIP(10,3)=  0.1000;   % CO2-F3


Cp_coeffs = zeros(nc, 4);
H_ig_ref  = zeros(nc, 1);

vt_method = 7;
Zc_c = 0.2905 - 0.085 * acentric;
vt_params = struct('Zc', Zc_c, 'components', ...
    {{'H2','N2','CO2','C1','C2','C3','C4C6','C7p1','C7p2','C7p3'}});

depth_range = [h_top, h_bot];


[GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = detect_sgoc( ...
    0, comp_ref, press_ref, temp_ref, h_ref, dTdh, depth_range, ...
    Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
    vt_method, vt_params);

%%
components = {'H2','N2','CO2','C1','C2','C3','C4C6','C7p1','C7p2','C7p3'};

[pressd, comp_liq] = pressdew_multicomp_ss(GOC_comp, 280.68e5, temp_ref, Pc, Tc, acentric, BIP, 10^(-8), 500, 'SRK', components);

fprintf('\nLiquid H2 at GOC: %.4f mol%%\n', comp_liq(1)*100);
fprintf('Vapour H2 at GOC: %.4f mol%%\n', GOC_comp(1)*100);

h_gas = linspace(h_top, GOC_depth, 80)';
ng    = length(h_gas);

xH2_gas = zeros(ng, 1);
P_gas   = zeros(ng, 1);
rho_gas = zeros(ng, 1);

for i = 1:ng
    [comp_h, P_h, T_h, ~, ~] = main_hasse(h_gas(i), h_ref, comp_ref, press_ref, temp_ref, ...
        dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    xH2_gas(i) = comp_h(1) * 100;
    P_gas(i)   = P_h / 1e5;
    try
        rho_gas(i) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, 9, vt_params);
    catch
        rho_gas(i) = NaN;
    end
end

h_liq_vec = linspace(GOC_depth, h_bot, 40)';
nl        = length(h_liq_vec);

xH2_liq  = zeros(nl, 1);
P_liq    = zeros(nl, 1);
rho_liq  = zeros(nl, 1);

for i = 1:nl
    [comp_l, P_l, T_l, ~, ~] = main_hasse(h_liq_vec(i), GOC_depth, comp_liq, ...
        GOC_pressure, GOC_temp, dTdh, Pc, Tc, acentric, BIP, M_gmol, ...
        Cp_coeffs, H_ig_ref, vt_method, vt_params);
    xH2_liq(i) = comp_l(1) * 100;
    P_liq(i)   = P_l / 1e5;
    try
        rho_liq(i) = calculate_density(comp_l, P_l, T_l, Pc, Tc, acentric, BIP, M_gmol, 9, vt_params);
    catch
        rho_liq(i) = NaN;
    end
end

fprintf('\nComputing P_dew profile along gas column...\n');
P_dew_gas = zeros(ng, 1);
comp_all  = zeros(ng, nc);

for i = 1:ng
    [comp_h, ~, ~, ~, ~] = main_hasse(h_gas(i), h_ref, comp_ref, press_ref, temp_ref, ...
        dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    comp_all(i,:) = comp_h';
end

for i = 1:ng
    try
        [psat_i, ~] = pressdew_multicomp_ss(comp_all(i,:)', P_gas(i)*1e5, temp_ref, ...
            Pc, Tc, acentric, BIP, 1e-8, 500, 'SRK', ...
            {'H2','N2','CO2','C1','C2','C3','C4C6','C7p1','C7p2','C7p3'});
        P_dew_gas(i) = psat_i / 1e5;
    catch
        P_dew_gas(i) = NaN;
    end
    if mod(i,10)==0
        fprintf('  Depth %.1f m: P_dew = %.2f bar\n', h_gas(i), P_dew_gas(i));
    end
end
fprintf('P_dew profile done.\n\n');

fprintf('Computing P_bub profile along liquid column...\n');
P_bub_liq = zeros(nl, 1);
comp_liq_all = zeros(nl, nc);

for i = 1:nl
    [comp_l, ~, ~, ~, ~] = main_hasse(h_liq_vec(i), GOC_depth, comp_liq, ...
        GOC_pressure, GOC_temp, dTdh, Pc, Tc, acentric, BIP, M_gmol, ...
        Cp_coeffs, H_ig_ref, vt_method, vt_params);
    comp_liq_all(i,:) = comp_l';
end

for i = 1:nl
    P_guesses = [P_liq(i)*1e5 * 1.1, P_liq(i)*1e5 * 1.5, P_liq(i)*1e5 * 2.0, GOC_pressure * 1.2];
    pbub_found = NaN;
    for ig = 1:length(P_guesses)
        try
            [pbub_i, ~] = pressbub_multicomp_ss(comp_liq_all(i,:)', P_guesses(ig), temp_ref, ...
                Pc, Tc, acentric, BIP, 1e-8, 500, 'SRK', ...
                {'H2','N2','CO2','C1','C2','C3','C4C6','C7p1','C7p2','C7p3'});
            if ~isnan(pbub_i) && pbub_i > 0
                pbub_found = pbub_i / 1e5;
                break;
            end
        catch
        end
    end
    P_bub_liq(i) = pbub_found;
    if mod(i,5)==0
        fprintf('  Depth %.1f m: P_bub = %.2f bar\n', h_liq_vec(i), P_bub_liq(i));
    end
end
fprintf('P_bub profile done. NaN count = %d / %d\n\n', sum(isnan(P_bub_liq)), nl);

%% ================================================================
%  Plot — three thesis-style panels
%% ================================================================
fn       = 'Times New Roman';
fs_tick  = 13;
fs_lab   = 14;
fs_pan   = 15;
fs_leg   = 12;
fs_zone  = 12;

pw    = 9;
ph    = 13;
gap_h = 1.5;
ml = 2.0; mr = 0.6; mt = 0.6; mb = 4.5;

col_vap   = [0.00 0.45 0.74];   % blue   - gas phase
col_liq   = [0.93 0.55 0.13];   % orange - incipient liquid
col_goc   = [0.20 0.20 0.20];   % near-black - GOC reference
col_shade = [0.95 0.85 0.70];   % light tan - two-phase zone

fig = figure;
fig.Units    = 'centimeters';
fig.Color    = 'w';
fig.Position = [2, 2, 3*pw + 2*gap_h + ml + mr, ph + mt + mb];
fig.PaperPositionMode = 'auto';

x_gas_data = {xH2_gas,  rho_gas,              P_gas};
xlabels    = {'H_{2} (mol%)', 'Density (kg/m^{3})', 'Pressure (bar)'};
panlabels  = {'(a)', '(b)', '(c)'};
x_liq_data = {xH2_liq, rho_liq, P_liq};

for k = 1:3
    left_cm = ml + (k-1)*(pw + gap_h);
    ax = axes('Units', 'centimeters', 'Position', [left_cm, mb, pw, ph]);
    hold on;

    all_x = [x_gas_data{k}; x_liq_data{k}];
    if k == 3
        all_x = [all_x; P_dew_gas; P_bub_liq(~isnan(P_bub_liq))];
    end
    x_lo = min(all_x);
    x_hi = max(all_x);
    x_pad = 0.04 * (x_hi - x_lo);
    if k == 1
        xl = [4.0, x_hi + x_pad];
    else
        xl = [x_lo - x_pad, x_hi + x_pad];
    end
    xlim(xl);

    hp = patch([xl(1) xl(2) xl(2) xl(1)], ...
        [GOC_depth GOC_depth h_bot h_bot], ...
        col_shade, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
    uistack(hp, 'bottom');

    h_vap = plot(x_gas_data{k}, h_gas, '-', ...
        'Color', col_vap, 'LineWidth', 2.6);

    h_liq = plot(x_liq_data{k}, h_liq_vec, '--', ...
        'Color', col_liq, 'LineWidth', 2.6);

    if k == 3
        h_dew = plot(P_dew_gas, h_gas, ':', ...
            'Color', col_vap, 'LineWidth', 1.8);
        h_bub = plot(P_bub_liq, h_liq_vec, ':', ...
            'Color', col_liq, 'LineWidth', 1.8);
        plot(GOC_pressure/1e5, GOC_depth, 'ko', ...
            'MarkerSize', 7, 'MarkerFaceColor', 'k');
    end

    plot(xl, [h_ref h_ref], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
    if k == 1
        text(xl(2) - 0.02*(xl(2)-xl(1)), h_ref - 1.5, ...
            sprintf('h_{ref} = %.1f m', h_ref), ...
            'FontName', fn, 'FontSize', fs_zone-1, 'Color', [0.4 0.4 0.4], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end

    plot(xl, [GOC_depth GOC_depth], '-', 'Color', col_goc, 'LineWidth', 1.0);
    if k == 1
        text(xl(2) - 0.02*(xl(2)-xl(1)), GOC_depth - 1.5, ...
            sprintf('GOC = %.1f m', GOC_depth), ...
            'FontName', fn, 'FontSize', fs_zone-1, 'Color', col_goc, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end

    if k == 1
        text(xl(1) + 0.04*(xl(2)-xl(1)), (GOC_depth+h_bot)/2, ...
            'Two-phase', 'FontName', fn, 'FontSize', fs_zone, ...
            'FontWeight', 'bold', ...
            'Color', [0.6 0.4 0.0], 'VerticalAlignment', 'middle');
    end

    hold off;
    set(ax, 'YDir', 'reverse', 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'FontName', fn, 'FontSize', fs_tick, 'Layer', 'top', ...
        'LineWidth', 0.9, 'TickLength', [0.015 0.015]);
    xlim(xl);

    xlabel(xlabels{k}, 'FontName', fn, 'FontSize', fs_lab);
    ylabel('Depth (m)',  'FontName', fn, 'FontSize', fs_lab);
    text(0.05, 0.97, panlabels{k}, 'Units', 'normalized', ...
        'FontName', fn, 'FontSize', fs_pan, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

end

%% Legend
ax_leg = axes('Units', 'centimeters', 'Position', [ml, 0.1, 3*pw+2*gap_h, mb-0.2], ...
    'Visible', 'off');
hold(ax_leg, 'on');

h1 = plot(ax_leg, NaN, NaN, '-',  'Color', col_vap, 'LineWidth', 2.6);
h2 = plot(ax_leg, NaN, NaN, '--', 'Color', col_liq, 'LineWidth', 2.6);
h3 = plot(ax_leg, NaN, NaN, ':',  'Color', col_vap, 'LineWidth', 1.8);
h4 = plot(ax_leg, NaN, NaN, ':',  'Color', col_liq, 'LineWidth', 1.8);
h5 = plot(ax_leg, NaN, NaN, '-',  'Color', col_goc, 'LineWidth', 1.0);

leg_labels = {
    'Gas phase profile', ...
    'Liquid phase profile', ...
    'Gas P_{dew}', ...
    'Liquid P_{bub}', ...
    sprintf('GOC at %.1f m', GOC_depth)
};

leg = legend(ax_leg, [h1, h2, h3, h4, h5], leg_labels, ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3, ...
    'FontName', fn, 'FontSize', fs_leg, ...
    'Box', 'on', 'EdgeColor', [0.4 0.4 0.4], ...
    'LineWidth', 0.7, 'Color', 'w');
leg.ItemTokenSize = [28, 14];
leg.Units = 'centimeters';
total_fig_w = 3*pw + 2*gap_h + ml + mr;
leg.Position(2) = 0.4;
leg.Position(1) = (total_fig_w - leg.Position(3)) / 2;


%% Volumetric H2 storage density calculation
M_H2 = M_gmol(1);

c_H2_gas = zeros(ng, 1);
for i = 1:ng
    M_mix_g = sum(comp_all(i,:)' .* M_gmol);
    c_H2_gas(i) = (comp_all(i,1) * M_H2 / M_mix_g) * rho_gas(i);
end

c_H2_liq = zeros(nl, 1);
for i = 1:nl
    M_mix_l = sum(comp_liq_all(i,:)' .* M_gmol);
    c_H2_liq(i) = (comp_liq_all(i,1) * M_H2 / M_mix_l) * rho_liq(i);
end

fprintf('\n=== Volumetric H2 Storage Density at GOC ===\n');
fprintf('Gas  phase: z_H2 = %.2f mol%%,  rho = %.1f kg/m3,  M_mix = %.1f g/mol,  c_H2 = %.3f kg H2/m3\n', ...
    comp_all(end,1)*100, rho_gas(end), sum(comp_all(end,:)'.*M_gmol), c_H2_gas(end));
fprintf('Liquid phase: z_H2 = %.2f mol%%,  rho = %.1f kg/m3,  M_mix = %.1f g/mol,  c_H2 = %.3f kg H2/m3\n', ...
    comp_liq_all(1,1)*100, rho_liq(1), sum(comp_liq_all(1,:)'.*M_gmol), c_H2_liq(1));
fprintf('Ratio gas/liquid: %.2f x\n', c_H2_gas(end) / c_H2_liq(1));