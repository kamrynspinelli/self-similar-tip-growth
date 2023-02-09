% flat tip, Guassian with \mu > 0
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight'); % set up the layout
nexttile;
load('secretion_to_shape_search_evolution_N_96_frames_700_flat.mat');
plot_splines(rv);
% xlim([(max(rv(1,:))+0.1 - 3.3*0.84*1.1) (max(rv(1,:))+0.1)]);
xlim([(max(rv(1,:))+0.1 - 2*0.84*1.1) (max(rv(1,:))+0.1)]);
ylim([0 0.84*1.1]);
pa = pbaspect;
nexttile;
hold on;
plot(0:0.001:s(end-1), gam(0:0.001:s(end-1)), 'LineWidth', 2.0);
plot(s(1:end-1), gam_num, 'LineWidth', 2.0);
load('secretion_to_shape_search_evolution_N_48_frames_700_flat.mat');
plot(s(1:end-1), gam_num, 'LineWidth', 2.0);
legend('Analytical \gamma: \mu/\sigma = 1.6', ['Reconstructed \gamma', newline, '(fine discretization)'], ['Reconstructed \gamma', newline, '(coarse discretization)'])
xlim([0 3.2]);
pbaspect(pa);
exportgraphics(gcf, 'media/secretion_to_shape_evolution_twopanel_flat.png');
close all;

% pointy tip, Gaussian with \mu = 0
close all;
t = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight'); % set up the layout
nexttile;
load('secretion_to_shape_search_evolution_N_96_frames_1000_pointy.mat');
plot_splines(rv);
xlim([(max(rv(1,:))+0.1 - 2*1.18*1.1) (max(rv(1,:))+0.1)]);
ylim([0 1.18*1.1]);
pa = pbaspect;
nexttile;
hold on;
plot(0:0.001:s(end-1), gam(0:0.001:s(end-1)), 'LineWidth', 2.0);
plot(s(1:end-1), gam_num, 'LineWidth', 2.0);
load('secretion_to_shape_search_evolution_N_48_frames_1000_pointy.mat');
plot(s(1:end-1), gam_num, 'LineWidth', 2.0);
legend('Analytical \gamma: \mu/\sigma = 0', ['Reconstructed \gamma', newline, '(fine discretization)'], ['Reconstructed \gamma', newline, '(coarse discretization)'])
xlim([0 3.2]);
ylim([-0.02 1]);
pbaspect(pa);
exportgraphics(gcf, 'media/secretion_to_shape_evolution_twopanel_pointy.png');
close all;