clc; clear; close all;

comp = [1.41; 0.140; 49.300; 9.780; 4.670; 0.620; 1.320; 0.220; 2.96; 1.82; ...
        3.05; 3.46; 2.69; 11.72; 4.51; 2.33];
z = comp ./ sum(comp);

Pc = [33.94; 73.76; 46; 48.84; 42.46; 36.48; 38; 33.84; 33.74; 29.69; ...
      31.95; 29.76; 26.67; 20.89; 15.71; 13.37] * 1e5;

Tc = [126.2; 304.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.4; ...
      535.334; 555.866; 577.173; 635.474; 736.020; 895.310];

acentric = [0.04; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; 0.296; ...
            0.4679; 0.4999; 0.5399; 0.6691; 0.9291; 1.2656];

ncomp = length(z);

BIP = zeros(ncomp, ncomp);
BIP(1,2) = -0.0315; BIP(2,1) = -0.0315;
BIP(1,3) =  0.0278; BIP(3,1) =  0.0278;
BIP(1,4) =  0.0407; BIP(4,1) =  0.0407;
BIP(1,5) =  0.0763; BIP(5,1) =  0.0763;
BIP(1,6) =  0.0944; BIP(6,1) =  0.0944;
BIP(1,7) =  0.07;   BIP(7,1) =  0.07;
BIP(1,8) =  0.0867; BIP(8,1) =  0.0867;
BIP(1,9) =  0.0878; BIP(9,1) =  0.0878;
for i = 10:ncomp
    BIP(1,i) = 0.08; BIP(i,1) = 0.08;
end
for i = 3:8
    BIP(2,i) = 0.12; BIP(i,2) = 0.12;
end
for i = 9:ncomp
    BIP(2,i) = 0.10; BIP(i,2) = 0.10;
end

comp_names = {'N2','CO2','C1','C2','C3','iC4','nC4','iC5','nC5','C6', ...
              'C7','C8','C9', 'C10-C15','C16-C25','C26-C36'};

T_crit = 315.09 + 273.15;
P_crit = 224.05e5;

%% Main envelope (bubble formulation traversing through CP)
env = pt_phase_envelope(z, Pc, Tc, acentric, BIP, 'SRK', comp_names, ...
    'plot',        false,     ...
    'max_steps',   1000,      ...
    'ds_init',     0.08,      ...
    'tol',         1e-9,      ...
    'maxiter',     40,        ...
    'T_range',     [250, T_crit+200], ...
    'formulation', 'bubble');

T_all = env.T_all - 273.15;
P_all = env.P_all / 1e5;

betas = [0.2, 0.4];
quality_curves = cell(length(betas), 1);

for bi = 1:length(betas)
    beta_val = betas(bi);
    fprintf('\n=== Computing quality line beta = %.1f ===\n', beta_val);
    env_q = pt_phase_envelope(z, Pc, Tc, acentric, BIP, 'SRK', comp_names, ...
        'plot',        false,     ...
        'max_steps',   1000,      ...
        'ds_init',     0.08,      ...
        'tol',         1e-9,      ...
        'maxiter',     40,        ...
        'T_range',     [250, T_crit+200], ...
        'formulation', 'quality', ...
        'beta',        beta_val);
    quality_curves{bi}.T = env_q.T_all - 273.15;
    quality_curves{bi}.P = env_q.P_all / 1e5;
end

%% Plotting
fig = figure('Units','centimeters','Position',[2 2 18 13],'Color','w', ...
    'PaperPositionMode','auto');
ax = axes('Parent', fig, 'Position',[0.11 0.13 0.66 0.80]); 
hold(ax,'on'); box(ax,'on');

q_colors = [0.93 0.55 0.13;     % orange  — V/F = 0.2
            0.18 0.69 0.34;     % green   — V/F = 0.4
            0.58 0.30 0.71;     % purple  — V/F = 0.6
            0.95 0.79 0.10];    % yellow  — V/F = 0.8

q_styles = {'--', '--', '--', '--'};
q_lw     = 1.6;

q_handles = gobjects(length(betas), 1);
for bi = 1:length(betas)
    q_handles(bi) = plot(ax, quality_curves{bi}.T, quality_curves{bi}.P, q_styles{bi}, ...
        'Color', q_colors(bi,:), 'LineWidth', q_lw, ...
        'DisplayName', sprintf('V/F = %.1f', betas(bi)));
end

h_env = plot(ax, T_all, P_all, '-', 'Color',[0.00 0.32 0.62], 'LineWidth', 2.6, ...
    'DisplayName','Phase envelope');

h_cp = plot(ax, T_crit-273.15, P_crit/1e5, 'p', ...
    'MarkerSize', 18, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.85 0.10 0.10], ...
    'LineWidth', 1.0, ...
    'DisplayName', sprintf('CP (%.1f\\circC, %.0f bar)', T_crit-273.15, P_crit/1e5));

set(ax, 'FontSize', 12, 'FontName','Times New Roman', ...
        'TickDir','in','XMinorTick','on','YMinorTick','on', ...
        'TickLength',[0.015 0.015],'LineWidth', 1.0);
xlabel(ax,'Temperature (\circC)','FontSize',13,'FontName','Times New Roman','FontWeight','bold');
ylabel(ax,'Pressure (bar)',       'FontSize',13,'FontName','Times New Roman','FontWeight','bold');

lg = legend(ax,[h_env; h_cp; q_handles], 'Location','eastoutside', ...
    'FontSize',10,'FontName','Times New Roman','Box','on', ...
    'EdgeColor',[0.4 0.4 0.4],'LineWidth',0.6);
lg.ItemTokenSize = [26, 12];

xlim(ax,[-50 450]); ylim(ax,[0 300]);

fprintf('\nEnvelope: %d points\n', length(T_all));
for bi = 1:length(betas)
    fprintf('V/F = %.1f: %d points\n', betas(bi), length(quality_curves{bi}.T));
end
fprintf('CP (input):    %.2f C, %.2f bar\n', T_crit-273.15, P_crit/1e5);