clear all;
a = 8;
stiffness = 24;

growth_anisotropic_hyphoid_extended_shank(a, stiffness);
growth_anisotropic_hyphoid(a, stiffness);

hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
load(['hyphoid_a_', num2str(a), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
plot(s, gam, 'LineWidth', 2.0); 
ylim([0 2.4]);
xlim([0 6]);
xlabel('s');
ylabel('\gamma');
ax = gca;
set(ax, 'fontsize', 12);
nexttile;
load(['hyphoid_extended_shank_a_', num2str(a), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); plot(s, gam, 'LineWidth', 2.0);
ylim([0 2.4]);
xlim([0 6]);
xlabel('s');
ylabel('\gamma');
ax = gca;
set(ax, 'fontsize', 12);

exportgraphics(gcf, ['media/hyphoid_extended_vs_normal_shank_secretion_a_', num2str(a), '_stiffness_', num2str(stiffness), '.png']);