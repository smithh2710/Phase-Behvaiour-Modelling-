clear; clc; close all;

R = 8.314462;
g = 9.80665;
h_ref     = 2600;
temp      = 352;
thickness = 122;
press_ref = 310e5;
dh        = 1;
tol       = 1e-12;
maxiter   = 300;

h_top  = h_ref - thickness/2;
h_bot  = h_ref + thickness/2;
h_down = (h_ref + dh : dh : h_bot)';
h_up   = (h_ref - dh : -dh : h_top)';

% =========================================================================
%  CASE 3 — H2 (10%) + Hassanpouryouzband 6-component NG
%  Components: H2 | N2 | CO2 | CH4 | C2H6 | C3-C5 pseudo
%  Original NG (mol%): N2=1.95, CO2=1.14, C1=83.60, C2=7.48
%                      C3=3.92, iC4=0.81, nC4=0.81, iC5=0.15, nC5=0.14
%  C3-C5 lump (sum=5.83 mol%) — mole-fraction-weighted averaging
% =========================================================================

c35_zi   = [3.92; 0.81; 0.81; 0.15; 0.14] / 5.83;
c35_M    = [44.097; 58.123; 58.123; 72.150; 72.150];
c35_Tc   = [369.83; 408.14; 425.12; 460.43; 469.70];
c35_Pc   = [42.48;  36.48;  37.96;  33.78;  33.75];
c35_w    = [0.1523; 0.1770; 0.2002; 0.2281; 0.2515];

M_C35  = sum(c35_zi .* c35_M);
Tc_C35 = sum(c35_zi .* c35_Tc);
Pc_C35 = sum(c35_zi .* c35_Pc);
w_C35  = sum(c35_zi .* c35_w);

% 10% H2 + 90% NG blend (NG sums to 100 mol%)
comp3 = [0.10; 1.95; 1.14; 83.60; 7.48; 5.83] .* [1; repmat(0.9/100, 5, 1)];
comp3(1) = 0.10;
comp3(2:end) = [1.95; 1.14; 83.60; 7.48; 5.83] * 0.9/100;

M3  = [2.016;  28.014; 44.010; 16.043; 30.070; M_C35 ];
Tc3 = [33.15;  126.20; 304.18; 190.60; 305.32; Tc_C35];
Pc3 = [12.90;  33.98;  73.80;  46.00;  48.72;  Pc_C35] * 1e5;
w3  = [-0.220; 0.037;  0.239;  0.011;  0.0995; w_C35 ];

BIP3 = zeros(6);
BIP3(1,2)=0.10;  BIP3(2,1)=0.10;
BIP3(1,3)=0.10;  BIP3(3,1)=0.10;
BIP3(1,4)=0.04;  BIP3(4,1)=0.04;
BIP3(1,5)=0.08;  BIP3(5,1)=0.08;
BIP3(1,6)=0.10;  BIP3(6,1)=0.10;
BIP3(2,3)=-0.017; BIP3(3,2)=-0.017;
BIP3(2,4)=0.031;  BIP3(4,2)=0.031;
BIP3(2,5)=0.052;  BIP3(5,2)=0.052;
BIP3(2,6)=0.085;  BIP3(6,2)=0.085;
BIP3(3,4)=0.120;  BIP3(4,3)=0.120;
BIP3(3,5)=0.130;  BIP3(5,3)=0.130;
BIP3(3,6)=0.135;  BIP3(6,3)=0.135;
BIP3(4,5)=0.0026; BIP3(5,4)=0.0026;
BIP3(4,6)=0.014;  BIP3(6,4)=0.014;
BIP3(5,6)=0.008;  BIP3(6,5)=0.008;

% =========================================================================
%  CASE 4 — H2 + Whitson gas condensate (11 components)
%  Components: H2 | C1 | C2 | C3 | iC4 | nC4 | iC5 | nC5 | C6 | C7-C9 | C10+
%
%  Whitson condensate (pre-split, mol%):
%    C1=70.24, C2=8.03, C3=4.40, iC4=0.94, nC4=1.88
%    iC5=0.84, nC5=0.90, C6=1.34, C7+=11.43 (MW=219, SG=0.839)
%
%  C7+ split (42% / 58% by mole, satisfies MW constraint 0.42*110+0.58*300=220≈219):
%    C7-C9 pseudo: z=4.80, MW=110, SG=0.765 → Tc=635K, Pc=21.3bar, ω=0.487
%    C10+  pseudo: z=6.63, MW=300, SG=0.885 → Tc=814K, Pc=13.3bar, ω=0.844
%  Critical properties computed via Soreide Tb + Kesler-Lee Tc/Pc + Edmister ω
% =========================================================================

wc_raw = [70.24; 8.03; 4.40; 0.94; 1.88; 0.84; 0.90; 1.34; 4.80; 6.63];
wc_ng  = wc_raw / sum(wc_raw);

comp4_10 = [0.10; wc_ng * 0.90];
comp4_20 = [0.20; wc_ng * 0.80];

M4  = [2.016;  16.043; 30.070; 44.097; 58.123; 58.123; 72.150; 72.150; 86.177; 110.0;  300.0 ];
Tc4 = [33.15;  190.60; 305.32; 369.83; 408.14; 425.12; 460.43; 469.70; 507.43; 635.0;  814.0 ];
Pc4 = [12.90;  46.00;  48.72;  42.48;  36.48;  37.96;  33.78;  33.75;  30.12;  21.3;   13.3  ] * 1e5;
w4  = [-0.220; 0.011;  0.0995; 0.1523; 0.1770; 0.2002; 0.2281; 0.2515; 0.3013; 0.487;  0.844 ];

nc4  = length(M4);
BIP4 = zeros(nc4);

h2_bip = [0, 0.04, 0.08, 0.10, 0.10, 0.11, 0.12, 0.12, 0.13, 0.14, 0.16];
c1_bip = [0, 0, 0.0026, 0.014, 0.016, 0.020, 0.022, 0.022, 0.026, 0.030, 0.038];
for j = 2:nc4
    BIP4(1,j) = h2_bip(j); BIP4(j,1) = h2_bip(j);
end
for j = 3:nc4
    BIP4(2,j) = c1_bip(j); BIP4(j,2) = c1_bip(j);
end
BIP4(3,4) = 0.005; BIP4(4,3) = 0.005;

% =========================================================================
%  STABILITY CHECK — flag if reference point is two-phase for Case 4
% =========================================================================
fprintf('=== STABILITY CHECK: Case 4, 10%% H2 at %.0f bar, %.0f K ===\n', ...
    press_ref/1e5, temp);
try
    [stab, ~] = stability_analysis_ssi(comp4_10, press_ref, temp, Pc4, Tc4, w4, BIP4, 'PR');
    if stab
        fprintf('  Single phase confirmed — grading calculation valid.\n\n');
    else
        fprintf('  WARNING: Two-phase at reference point — single-phase model not valid here.\n\n');
    end
catch
    fprintf('  stability_analysis_ssi unavailable — check manually via flash at ref point.\n\n');
end

% =========================================================================
%  RUN ISOTHERMAL GRADING
% =========================================================================
fprintf('Running Case 3 ...\n');
[h3, z3, ~, ~] = run_isothermal_grading(comp3, press_ref, temp, ...
    Pc3, Tc3, w3, M3*1e-3, BIP3, h_ref, h_down, h_up, dh, R, g, tol, maxiter);

fprintf('Running Case 4 (10%% H2) ...\n');
[h4_10, z4_10, ~, ~] = run_isothermal_grading(comp4_10, press_ref, temp, ...
    Pc4, Tc4, w4, M4*1e-3, BIP4, h_ref, h_down, h_up, dh, R, g, tol, maxiter);

fprintf('Running Case 4 (20%% H2) ...\n');
[h4_20, z4_20, ~, ~] = run_isothermal_grading(comp4_20, press_ref, temp, ...
    Pc4, Tc4, w4, M4*1e-3, BIP4, h_ref, h_down, h_up, dh, R, g, tol, maxiter);

fprintf('\n--- Delta H2 (top minus bottom) ---\n');
fprintf('  Case 3, 10%% H2 + lean NG:         %.5f mol%%\n', (z3(1,1)    - z3(end,1))   *100);
fprintf('  Case 4, 10%% H2 + condensate:      %.5f mol%%\n', (z4_10(1,1) - z4_10(end,1))*100);
fprintf('  Case 4, 20%% H2 + condensate:      %.5f mol%%\n', (z4_20(1,1) - z4_20(end,1))*100);

% =========================================================================
%  FIGURE 1 — H2 mol% vs depth (all three runs)
% =========================================================================
fn = 'Times New Roman'; fs = 14; lw = 1.8;

fig1 = figure;
fig1.Units = 'centimeters'; fig1.Color = 'w';
fig1.Position = [2 2 14 16];

ax1 = axes('Units','normalized','Position',[0.13 0.17 0.82 0.78]);
hold on;
hL = gobjects(3,1);
hL(1) = plot(z3(:,1)*100,    h3,    'Color',[0.00 0.45 0.70],'LineStyle','-', 'LineWidth',lw);
hL(2) = plot(z4_10(:,1)*100, h4_10, 'Color',[0.80 0.00 0.00],'LineStyle','--','LineWidth',lw);
hL(3) = plot(z4_20(:,1)*100, h4_20, 'Color',[0.90 0.60 0.00],'LineStyle','-.','LineWidth',lw);
hold off;
set(ax1,'YDir','reverse','Box','on','TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','FontName',fn,'FontSize',fs,'Layer','top');
xlabel('H_{2} (mol%)', 'FontName',fn,'FontSize',fs);
ylabel('Depth (m)',    'FontName',fn,'FontSize',fs);
text(0.05,0.97,'(a)','Units','normalized','FontName',fn,'FontSize',fs, ...
    'FontWeight','bold','VerticalAlignment','top');
leg1 = legend(hL, {'10% H_2 + lean NG', '10% H_2 + condensate', '20% H_2 + condensate'}, ...
    'Orientation','vertical','FontName',fn,'FontSize',fs-2,'Box','off','Location','best');

% =========================================================================
%  FIGURE 2 — All component profiles, Case 4 (10% H2 + Whitson condensate)
% =========================================================================
cnames = {'H_2','C_1','C_2','C_3','iC_4','nC_4','iC_5','nC_5','C_6','C_7-C_9','C_{10+}'};
ccolors = { [0.00 0.45 0.70], [0.00 0.60 0.50], [0.90 0.60 0.00], [0.80 0.00 0.00], ...
            [0.50 0.00 0.50], [0.35 0.70 0.90], [0.60 0.40 0.20], [0.40 0.40 0.40], ...
            [0.10 0.50 0.10], [0.95 0.45 0.00], [0.20 0.00 0.70] };
clstyles = {'-','--','-.', ':' ,'-','--','-.', ':' ,'-','--','-.'};

fig2 = figure;
fig2.Units = 'centimeters'; fig2.Color = 'w';
fig2.Position = [18 2 14 16];

ax2 = axes('Units','normalized','Position',[0.13 0.17 0.82 0.78]);
hold on;
hL2 = gobjects(nc4,1);
for i = 1:nc4
    hL2(i) = plot(z4_10(:,i)*100, h4_10, ...
        'Color',ccolors{i},'LineStyle',clstyles{i},'LineWidth',1.5);
end
hold off;
set(ax2,'YDir','reverse','Box','on','TickDir','in', ...
    'XMinorTick','on','YMinorTick','on','FontName',fn,'FontSize',fs,'Layer','top');
xlabel('Component (mol%)', 'FontName',fn,'FontSize',fs);
ylabel('Depth (m)',        'FontName',fn,'FontSize',fs);
text(0.05,0.97,'(b)','Units','normalized','FontName',fn,'FontSize',fs, ...
    'FontWeight','bold','VerticalAlignment','top');
legend(hL2, cnames, 'Orientation','vertical','FontName',fn,'FontSize',fs-3, ...
    'Box','off','Location','best');

% =========================================================================
%  LOCAL FUNCTION
% =========================================================================
function [all_h, all_z, all_P, all_rho] = run_isothermal_grading( ...
    comp_ref, press_ref, temp, Pc, Tc, w, M_kgmol, BIP, ...
    h_ref, h_down, h_up, dh, R, g, tol, maxiter)

    nc     = length(comp_ref);
    n_down = length(h_down);
    n_up   = length(h_up);

    comp_down = zeros(n_down, nc); P_down = zeros(n_down,1); rho_down = zeros(n_down,1);
    comp_up   = zeros(n_up,   nc); P_up   = zeros(n_up,  1); rho_up   = zeros(n_up,  1);

    [phi_ref, Z_ref] = fugacitycoef_multicomp_vapor(comp_ref, press_ref, temp, Pc, Tc, w, BIP, 'PR');
    rho_ref = press_ref * sum(comp_ref .* M_kgmol) / (Z_ref * R * temp);

    z_prev = comp_ref; P_prev = press_ref; phi_prev = phi_ref;
    for k = 1:n_down
        [~, Z_prev] = fugacitycoef_multicomp_vapor(z_prev, P_prev, temp, Pc, Tc, w, BIP, 'PR');
        rho_prev  = P_prev * sum(z_prev .* M_kgmol) / (Z_prev * R * temp);
        P_guess   = P_prev + rho_prev * g * dh;
        fi_target = (phi_prev .* z_prev * P_prev) .* exp(M_kgmol * g * dh / (R * temp));
        z_new = z_prev; P_new = P_guess;
        for iter = 1:maxiter
            [phi_new, Z_new] = fugacitycoef_multicomp_vapor(z_new, P_new, temp, Pc, Tc, w, BIP, 'PR');
            z_calc = fi_target ./ (phi_new * P_new);
            z_calc = z_calc / sum(z_calc);
            rho_new = P_new * sum(z_calc .* M_kgmol) / (Z_new * R * temp);
            P_calc  = P_prev + rho_new * g * dh;
            if max(abs(z_calc - z_new)) < tol && abs(P_calc - P_new)/P_new < tol; break; end
            z_new = z_calc; P_new = P_calc;
        end
        comp_down(k,:) = z_new'; P_down(k) = P_new;
        rho_down(k)    = P_new * sum(z_new .* M_kgmol) / (Z_new * R * temp);
        z_prev = z_new; P_prev = P_new; phi_prev = phi_new;
    end

    z_prev = comp_ref; P_prev = press_ref;
    [phi_prev, ~] = fugacitycoef_multicomp_vapor(comp_ref, press_ref, temp, Pc, Tc, w, BIP, 'PR');
    for k = 1:n_up
        [~, Z_prev] = fugacitycoef_multicomp_vapor(z_prev, P_prev, temp, Pc, Tc, w, BIP, 'PR');
        rho_prev  = P_prev * sum(z_prev .* M_kgmol) / (Z_prev * R * temp);
        P_guess   = P_prev - rho_prev * g * dh;
        fi_target = (phi_prev .* z_prev * P_prev) .* exp(-M_kgmol * g * dh / (R * temp));
        z_new = z_prev; P_new = P_guess;
        for iter = 1:maxiter
            [phi_new, Z_new] = fugacitycoef_multicomp_vapor(z_new, P_new, temp, Pc, Tc, w, BIP, 'PR');
            z_calc = fi_target ./ (phi_new * P_new);
            z_calc = z_calc / sum(z_calc);
            rho_new = P_new * sum(z_calc .* M_kgmol) / (Z_new * R * temp);
            P_calc  = P_prev - rho_new * g * dh;
            if max(abs(z_calc - z_new)) < tol && abs(P_calc - P_new)/P_new < tol; break; end
            z_new = z_calc; P_new = P_calc;
        end
        comp_up(k,:) = z_new'; P_up(k) = P_new;
        rho_up(k)    = P_new * sum(z_new .* M_kgmol) / (Z_new * R * temp);
        z_prev = z_new; P_prev = P_new; phi_prev = phi_new;
    end

    all_h   = [flipud(h_up); h_ref; h_down];
    all_z   = [flipud(comp_up); comp_ref'; comp_down];
    all_P   = [flipud(P_up);  press_ref; P_down];
    all_rho = [flipud(rho_up); rho_ref;  rho_down];
end