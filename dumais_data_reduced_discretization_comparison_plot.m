clear all;
profile_num = 3;
stiffness = 24;

% growth_anisotropic_hyphoid_extended_shank(a, stiffness);
% growth_anisotropic_hyphoid(a, stiffness);

hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
load(['../cell-profiles/dumais-root-hair-', num2str(profile_num), '.mat']); 
plot(s, gam, 'LineWidth', 2.0); 
ylim([0 1.7]);
xlim([0 4]);
xlabel('s');
ylabel('\gamma');
ax = gca;
set(ax, 'fontsize', 12);
nexttile;
load(['../cell-profiles/dumais-root-hair-reduced-discretization', num2str(profile_num), '.mat']); plot(s, gam, 'LineWidth', 2.0);
ylim([0 1.7]);
xlim([0 4]);
xlabel('s');
ylabel('\gamma');
ax = gca;
set(ax, 'fontsize', 12);

exportgraphics(gcf, ['media/dumais_data_normal_vs_reduced_discretization_profile_', num2str(profile_num), '_stiffness_', num2str(stiffness), '.png']);