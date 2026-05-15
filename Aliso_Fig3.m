clc; clear; close all;

R = 8.3144598;

h_top    = 2290.3;
h_bot    = 2412.3;
h_ref    = 2351.3;
temp_ref = 352.32;
dTdh     = 0;
tol      = 1e-6;
maxiter  = 500;

h_vec = linspace(h_top, h_bot, 100)';
nd    = length(h_vec);

cases = struct();

% =========================================================================
%  Cases 1-2: multicomponent — Pc in Pa, press in Pa
%  main_hasse returns P_h in Pa → divide by 1e5 for bar
% =========================================================================
cases(1).label    = 'CH_4 + C_2H_6 cushion';
cases(1).comp     = [0.10; 0.85; 0.05];
cases(1).M        = [2.016; 16.043; 30.070];
cases(1).Tc       = [33.15; 190.60; 305.32];
cases(1).Pc       = [12.90; 46.00;  48.72] * 1e5;
cases(1).w        = [-0.220; 0.011; 0.0995];
cases(1).BIP      = zeros(3,3);
cases(1).cnames   = {'H2','CH4','C2H6'};
cases(1).press    = 310.26e5;
cases(1).Pa_units = true;

% DZT depleted vapor (Zhang et al. 2020) + 10% H2
% Pseudo-component critical properties via Twu (1984) / Lee-Kesler
% C4-C6: MW=66.9, C7+1: MW=107.8, C7+2: MW=198.5, C7+3: MW=335.1
dzt_dep = [1.098; 1.917; 63.477; 7.821; 5.049; 7.443; 2.781; 0.396; 0.018] / 100;
comp2   = [0.10; dzt_dep];

M2  = [2.016;  44.010; 28.014; 16.043; 30.070; 44.097; 66.900; 107.800; 198.500; 335.100];
Tc2 = [33.20;  304.18; 126.20; 190.60; 305.32; 369.83; 460.00; 578.000; 718.000; 868.000];
Pc2 = [12.97;  73.76;  33.94;  46.00;  48.84;  42.46;  33.50;  24.800;  16.200;  11.800] * 1e5;
w2  = [-0.220; 0.225;  0.040;  0.008;  0.098;  0.152;  0.235;  0.412;   0.723;   1.050];

nc2  = 10;
BIP2 = zeros(nc2);

% Component index: 1=H2, 2=CO2, 3=N2, 4=C1, 5=C2, 6=C3, 7=iC4-nC6, 8=F1, 9=F2, 10=F3
% Values read directly from PVTsim Interact Param table (Zhang et al. 2020 fluid)

% H2 BIPs
BIP2(1,2) = -0.3426; BIP2(2,1) = -0.3426;  % H2-CO2
BIP2(1,3) =  0.0069; BIP2(3,1) =  0.0069;  % H2-N2
BIP2(1,4) = -0.0220; BIP2(4,1) = -0.0220;  % H2-C1
BIP2(1,5) = -0.1667; BIP2(5,1) = -0.1667;  % H2-C2
BIP2(1,6) = -0.2359; BIP2(6,1) = -0.2359;  % H2-C3
BIP2(1,7) = -0.5071; BIP2(7,1) = -0.5071;  % H2-iC4nC6
BIP2(1,8) =  0.2500; BIP2(8,1) =  0.2500;  % H2-F1
BIP2(1,9) =  0.2500; BIP2(9,1) =  0.2500;  % H2-F2
BIP2(1,10)=  0.2500; BIP2(10,1)=  0.2500;  % H2-F3

% CO2 BIPs
BIP2(2,3) = -0.0315; BIP2(3,2) = -0.0315;  % CO2-N2
BIP2(2,4) =  0.1200; BIP2(4,2) =  0.1200;  % CO2-C1
BIP2(2,5) =  0.1200; BIP2(5,2) =  0.1200;  % CO2-C2
BIP2(2,6) =  0.1200; BIP2(6,2) =  0.1200;  % CO2-C3
BIP2(2,7) =  0.1200; BIP2(7,2) =  0.1200;  % CO2-iC4nC6
BIP2(2,8) =  0.1000; BIP2(8,2) =  0.1000;  % CO2-F1
BIP2(2,9) =  0.1000; BIP2(9,2) =  0.1000;  % CO2-F2
BIP2(2,10)=  0.1000; BIP2(10,2)=  0.1000;  % CO2-F3

% N2 BIPs
BIP2(3,4) =  0.0278; BIP2(4,3) =  0.0278;  % N2-C1
BIP2(3,5) =  0.0407; BIP2(5,3) =  0.0407;  % N2-C2
BIP2(3,6) =  0.0763; BIP2(6,3) =  0.0763;  % N2-C3
BIP2(3,7) =  0.0702; BIP2(7,3) =  0.0702;  % N2-iC4nC6
BIP2(3,8) =  0.0800; BIP2(8,3) =  0.0800;  % N2-F1
BIP2(3,9) =  0.0800; BIP2(9,3) =  0.0800;  % N2-F2
BIP2(3,10)=  0.0800; BIP2(10,3)=  0.0800;  % N2-F3

% All remaining pairs (C1-C2, C1-C3 ... F1-F2 etc.) = 0.0 (default zeros)

cases(2).label    = 'Residual(depleted gas condensate)';
cases(2).comp     = comp2;  cases(2).M  = M2;   cases(2).Tc = Tc2;
cases(2).Pc       = Pc2;    cases(2).w  = w2;   cases(2).BIP = BIP2;
cases(2).cnames   = {'H2','CO2','N2','C1','C2','C3','C4C6','C7p1','C7p2','C7p3'};
cases(2).press    = 310.26e5;
cases(2).Pa_units = true;

% =========================================================================
%  Cases 3-6: pure binaries — Pc in bar, press in bar
%  main_hasse returns P_h in bar → use directly
% =========================================================================
cases(3).label    = 'CH_4 cushion';
cases(3).comp     = [0.10; 0.90];
cases(3).M        = [2.016; 16.043];
cases(3).Tc       = [33.20; 190.60];
cases(3).Pc       = [12.97; 46.00];
cases(3).w        = [-0.220; 0.008];
BIP3 = zeros(2); BIP3(1,2) = -0.0220; BIP3(2,1) = -0.0220;
cases(3).BIP      = BIP3;
cases(3).cnames   = {'H2','CH4'};
cases(3).press    = 310.26;
cases(3).Pa_units = false;

cases(4).label    = 'C_2H_6 cushion';
cases(4).comp     = [0.10; 0.90];
cases(4).M        = [2.016; 30.070];
cases(4).Tc       = [33.20; 305.32];
cases(4).Pc       = [12.97; 48.84];
cases(4).w        = [-0.220; 0.0980];
BIP4 = zeros(2); BIP4(1,2) = -0.1667; BIP4(2,1) = -0.1667;
cases(4).BIP      = BIP4;
cases(4).cnames   = {'H2','C2H6'};
cases(4).press    = 310.26;
cases(4).Pa_units = false;

cases(5).label    = 'Pure N_2 cushion';
cases(5).comp     = [0.10; 0.90];
cases(5).M        = [2.016; 28.014];
cases(5).Tc       = [33.20; 126.20];
cases(5).Pc       = [12.97; 33.94];
cases(5).w        = [-0.220; 0.047];
BIP5 = zeros(2); BIP5(1,2) = 0.0069; BIP5(2,1) = 0.0069;
cases(5).BIP      = BIP5;
cases(5).cnames   = {'H2','N2'};
cases(5).press    = 310.26;
cases(5).Pa_units = false;

cases(6).label    = 'CO_2 cushion';
cases(6).comp     = [0.10; 0.90];
cases(6).M        = [2.016; 44.010];
cases(6).Tc       = [33.20; 304.18];
cases(6).Pc       = [12.97; 73.76];
cases(6).w        = [-0.220; 0.225];
BIP6 = zeros(2); BIP6(1,2) = -0.3426; BIP6(2,1) = -0.3426;
cases(6).BIP      = BIP6;
cases(6).cnames   = {'H2','CO2'};
cases(6).press    = 310.26;
cases(6).Pa_units = false;

ncases = length(cases);

xH2_all = zeros(nd, ncases);
P_all   = zeros(nd, ncases);
rho_all = zeros(nd, ncases);

for c = 1:ncases
    nc_c  = length(cases(c).comp);
    Zc_c  = 0.2905 - 0.085 * cases(c).w;
    Cp_c  = zeros(nc_c, 4);
    Hig_c = zeros(nc_c, 1);
    vtp   = struct('Zc', Zc_c, 'components', {cases(c).cnames});

    for i = 1:nd
        [comp_h, P_h, T_h, ~, ~] = main_hasse(h_vec(i), h_ref, cases(c).comp, ...
            cases(c).press, temp_ref, dTdh, cases(c).Pc, cases(c).Tc, cases(c).w, ...
            cases(c).BIP, cases(c).M, Cp_c, Hig_c, 7, vtp);
        xH2_all(i,c) = comp_h(1) * 100;
        if cases(c).Pa_units
            P_all(i,c) = P_h / 1e5;
        else
            P_all(i,c) = P_h;
        end
        try
            rho_all(i,c) = calculate_density(comp_h, P_h, T_h, cases(c).Pc, cases(c).Tc, ...
                cases(c).w, cases(c).BIP, cases(c).M, 9, vtp);
        catch
            rho_all(i,c) = NaN;
        end
    end
    fprintf('Case %d (%s): xH2 top=%.4f  bot=%.4f  Delta=%.4f mol%%\n', ...
        c, cases(c).label, xH2_all(1,c), xH2_all(end,c), xH2_all(1,c)-xH2_all(end,c));
end

ref_idx   = round(nd / 2);
rho_ratio = rho_all ./ rho_all(ref_idx, :);

fn       = 'Times New Roman';
fs_tick  = 13;
fs_lab   = 14;
fs_pan   = 15;
fs_leg   = 12;

pw    = 9;
ph    = 13;
gap_h = 1.5;
ml = 2.0; mr = 0.6; mt = 0.6; mb = 4.5;

line_colors = {
    [0.85 0.10 0.10],   % red         - Aliso Canyon (operational benchmark)
    [0.00 0.45 0.74],   % blue        - DZT depleted gas condensate (novel result)
    [0.00 0.60 0.50],   % teal        - Pure CH4
    [0.93 0.69 0.13],   % amber       - Pure C2H6
    [0.58 0.30 0.71],   % purple      - Pure N2
    [0.93 0.55 0.13]    % orange      - Pure CO2
};
line_styles = {'-', '--', '-.', ':', '--', '-.'};
line_widths = [2.6, 2.6, 2.0, 2.0, 2.0, 2.0];

fig = figure;
fig.Units    = 'centimeters';
fig.Color    = 'w';
fig.Position = [2, 2, 3*pw + 2*gap_h + ml + mr, ph + mt + mb];
fig.PaperPositionMode = 'auto';

x_data    = {xH2_all, rho_ratio, P_all};
xlabels   = {'H_{2} (mol%)', '\rho / \rho_{ref} (-)', 'Pressure (bar)'};
panlabels = {'(a)', '(b)', '(c)'};

for k = 1:3
    left_cm = ml + (k-1)*(pw + gap_h);
    ax = axes('Units', 'centimeters', 'Position', [left_cm, mb, pw, ph]);
    hold on;
    h_lines = gobjects(ncases, 1);
    for c = 1:ncases
        h_lines(c) = plot(x_data{k}(:,c), h_vec, ...
            'LineStyle', line_styles{c}, ...
            'Color',     line_colors{c}, ...
            'LineWidth', line_widths(c));
    end
    yline(h_ref, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    hold off;
    set(ax, 'YDir', 'reverse', 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'FontName', fn, 'FontSize', fs_tick, 'Layer', 'top', ...
        'LineWidth', 0.9, 'TickLength', [0.015 0.015]);
    xlabel(xlabels{k}, 'FontName', fn, 'FontSize', fs_lab);
    ylabel('Depth (m)', 'FontName', fn, 'FontSize', fs_lab);
    text(0.05, 0.97, panlabels{k}, 'Units', 'normalized', ...
        'FontName', fn, 'FontSize', fs_pan, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    if k == 2
        case_labels = {cases.label};
        leg = legend(h_lines, case_labels, ...
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
    end
end