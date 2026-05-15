clc; clear; close all;

comp = [1.41; 0.140; 49.300; 9.780; 4.670; 0.620; 1.320; 0.220; 2.96; 1.82; ...
        3.05; 3.46; 2.69; 11.72; 4.51; 2.33];
comp = comp ./ sum(comp);

Pc = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; ...
      31.95; 29.76; 26.67; 20.89; 15.71; 13.37] * 1e5;

Tc = [126.2; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.4; ...
      535.334; 555.866; 577.173; 635.474; 736.020; 895.310];

acentric = [0.04; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; 0.296; ...
            0.4679; 0.4999; 0.5399; 0.6691; 0.9291; 1.2656];

M_gmol = [28.014; 44.01; 16.043; 30.07; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; ...
          96; 107; 121; 165.41; 269.829; 477.892];

ncomp = length(comp);

BIP = zeros(ncomp, ncomp);
BIP(1,2) = -0.0315; BIP(2,1) = -0.0315;
BIP(1,3) = 0.0278;  BIP(3,1) = 0.0278;
BIP(1,4) = 0.0407;  BIP(4,1) = 0.0407;
BIP(1,5) = 0.0763;  BIP(5,1) = 0.0763;
BIP(1,6) = 0.0944;  BIP(6,1) = 0.0944;
BIP(1,7) = 0.07;    BIP(7,1) = 0.07;
BIP(1,8) = 0.0867;  BIP(8,1) = 0.0867;
BIP(1,9) = 0.0878;  BIP(9,1) = 0.0878;
for i = 10:ncomp
    BIP(1,i) = 0.08; BIP(i,1) = 0.08;
end
for i = 3:8
    BIP(2,i) = 0.12; BIP(i,2) = 0.12;
end
for i = 9:ncomp
    BIP(2,i) = 0.10; BIP(i,2) = 0.10;
end

Zc = 0.2905 - 0.085 .* acentric;
components = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6', ...
              'C7','C8','C9','C10-C15','C16-C25','C26-C36'};

h_ref     = 2000;
dTdh      = 0;
vt_method = 9;
vt_params = struct('Zc', Zc, 'components', {components});

Cp_coeffs = zeros(ncomp, 4);
H_ig_ref  = zeros(ncomp, 1);

T_refs = [65; 70; 75; 80; 85] + 273.15;
P_refs = [265; 285; 305; 325; 345] * 1e5;
labels = {'T = 65°C, P = 265 bar', ...
          'T = 70°C, P = 285 bar', ...
          'T = 75°C, P = 305 bar', ...
          'T = 80°C, P = 325 bar', ...
          'T = 85°C, P = 345 bar'};

depths  = (2000:-50:1550)';
nd      = length(depths);
nref    = 5;
idx_C1  = 3;
idx_C7p = 11:ncomp;

colors = [0.000 0.447 0.741;
          0.850 0.325 0.098;
          0.466 0.674 0.188;
          0.494 0.184 0.556;
          0.929 0.694 0.125];
lstyles = {'-', '--', '-.', ':', '-'};
lw = 1.5;

C1_all  = NaN(nd, nref);
C7p_all = NaN(nd, nref);
rho_all = NaN(nd, nref);

for r = 1:nref
    T_ref = T_refs(r);
    P_ref = P_refs(r);
    for i = 1:nd
        try
            [comp_h, P_h, T_h, ~, ~] = main_hasse(depths(i), h_ref, comp, P_ref, T_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
            C1_all(i,r)  = comp_h(idx_C1) * 100;
            C7p_all(i,r) = sum(comp_h(idx_C7p)) * 100;
            try
                rho_all(i,r) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
            catch
            end
        catch
        end
    end
end

rho_ratio = NaN(nd, nref);
for r = 1:nref
    if ~isnan(rho_all(1,r)) && rho_all(1,r) > 0
        rho_ratio(:,r) = rho_all(:,r) / rho_all(1,r);
    end
end

fmt_ax = @(ax) set(ax, ...
    'YDir',       'reverse', ...
    'FontSize',   13, ...
    'FontName',   'Times New Roman', ...
    'TickDir',    'in', ...
    'Box',        'on', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', [0.015 0.015], ...
    'LineWidth',  0.9);

fig_w = 14;
fig_h = 10;
lg_fs = 12;

%% Figure 1 — C1 mol%
fig1 = figure('Units','centimeters','Position',[2 14 fig_w fig_h],'Color','w');
ax1  = axes('Parent', fig1);
hold(ax1,'on');
for r = 1:nref
    plot(ax1, C1_all(:,r), depths, lstyles{r}, ...
        'Color', colors(r,:), 'LineWidth', lw, 'DisplayName', labels{r});
end
fmt_ax(ax1);
xlabel(ax1, 'C_1 (mol%)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax1, 'Depth (m)',  'FontSize', 14, 'FontName', 'Times New Roman');;
lg1 = legend(ax1, 'Location', 'southwest', 'FontSize', lg_fs, ...
    'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', [0.4 0.4 0.4], ...
    'LineWidth', 0.7, 'Color', 'w');
lg1.ItemTokenSize = [28, 14];
set(fig1, 'PaperPositionMode', 'auto');

%% Figure 2 — C7+ mol%
fig2 = figure('Units','centimeters','Position',[2 2 fig_w fig_h],'Color','w');
ax2  = axes('Parent', fig2);
hold(ax2,'on');
for r = 1:nref
    plot(ax2, C7p_all(:,r), depths, lstyles{r}, ...
        'Color', colors(r,:), 'LineWidth', lw, 'DisplayName', labels{r});
end
fmt_ax(ax2);
xlabel(ax2, 'C_{7+} (mol%)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax2, 'Depth (m)',     'FontSize', 14, 'FontName', 'Times New Roman');
lg2 = legend(ax2, 'Location', 'northwest', 'FontSize', lg_fs, ...
    'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', [0.4 0.4 0.4], ...
    'LineWidth', 0.7, 'Color', 'w');
lg2.ItemTokenSize = [28, 14];
set(fig2, 'PaperPositionMode', 'auto');

%% Figure 3 — Density ratio
fig3 = figure('Units','centimeters','Position',[18 2 fig_w fig_h],'Color','w');
ax3  = axes('Parent', fig3);
hold(ax3,'on');
for r = 1:nref
    v = ~isnan(rho_ratio(:,r));
    plot(ax3, rho_ratio(v,r), depths(v), lstyles{r}, ...
        'Color', colors(r,:), 'LineWidth', lw, 'DisplayName', labels{r});
end
xline(ax3, 1.0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
fmt_ax(ax3);
xlabel(ax3, '\rho / \rho_{ref}  (–)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel(ax3, 'Depth (m)',              'FontSize', 14, 'FontName', 'Times New Roman');
lg3 = legend(ax3, 'Location', 'southwest', 'FontSize', lg_fs, ...
    'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', [0.4 0.4 0.4], ...
    'LineWidth', 0.7, 'Color', 'w');
lg3.ItemTokenSize = [28, 14];
set(fig3, 'PaperPositionMode', 'auto');