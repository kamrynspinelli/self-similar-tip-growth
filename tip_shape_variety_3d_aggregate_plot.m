%% load the ellipsoid and hyphoid trendlines
% load('tip_variety_data_k_6_mu_6.mat');
load('tip_variety_data_k_6_mu_6_p1.mat');
a_nl_6 = a;
b_nl_6 = b;
gam_peaks_nl_6 = gam_peaks;
R0_nl_6 = R0;
RN_nl_6 = RN;
tip_lengths_nl_6 = tip_lengths;
ZN_nl_6 = ZN;
W_ints_nl_6 = W_ints;
W_nl_6 = W;
gams_nl_6 = gams;
L0Splines_nl_6 = L0Splines;
strainls_nl_6 = strainls;
load('tip_variety_data_k_6_mu_6_p2.mat');
a_nl_6 = [a_nl_6 a];
b_nl_6 = [b_nl_6 b];
gam_peaks_nl_6 = [gam_peaks_nl_6 gam_peaks];
R0_nl_6 = [R0_nl_6 R0];
RN_nl_6 = [RN_nl_6 RN];
tip_lengths_nl_6 = [tip_lengths_nl_6 tip_lengths];
ZN_nl_6 = [ZN_nl_6 ZN];
W_ints_nl_6 = [W_ints_nl_6; W_ints];
W_nl_6 = [W_nl_6 W];
gams_nl_6 = [gams_nl_6; gams];
L0Splines_nl_6 = [L0Splines_nl_6; L0Splines];
strainls_nl_6 = [strainls_nl_6; strainls];

% load('tip_variety_data_k_12_mu_12.mat');
load('tip_variety_data_k_12_mu_12_p1.mat');
a_nl_12 = a;
b_nl_12 = b;
gam_peaks_nl_12 = gam_peaks;
R0_nl_12 = R0;
RN_nl_12 = RN;
tip_lengths_nl_12 = tip_lengths;
ZN_nl_12 = ZN;
W_ints_nl_12 = W_ints;
W_nl_12 = W;
gams_nl_12 = gams;
L0Splines_nl_12 = L0Splines;
strainls_nl_12 = strainls;
load('tip_variety_data_k_12_mu_12_p2.mat');
a_nl_12 = [a_nl_12 a];
b_nl_12 = [b_nl_12 b];
gam_peaks_nl_12 = [gam_peaks_nl_12 gam_peaks];
R0_nl_12 = [R0_nl_12 R0];
RN_nl_12 = [RN_nl_12 RN];
tip_lengths_nl_12 = [tip_lengths_nl_12 tip_lengths];
ZN_nl_12 = [ZN_nl_12 ZN];
W_ints_nl_12 = [W_ints_nl_12; W_ints];
W_nl_12 = [W_nl_12 W];
gams_nl_12 = [gams_nl_12; gams];
L0Splines_nl_12 = [L0Splines_nl_12; L0Splines];
strainls_nl_12 = [strainls_nl_12; strainls];

load('tip_variety_data_k_24_mu_24.mat');
a_nl_24 = a;
b_nl_24 = b;
gam_peaks_nl_24 = gam_peaks;
R0_nl_24 = R0;
RN_nl_24 = RN;
tip_lengths_nl_24 = tip_lengths;
ZN_nl_24 = ZN;
W_ints_nl_24 = W_ints;
W_nl_24 = W;
gams_nl_24 = gams;
L0Splines_nl_24 = L0Splines;
strainls_nl_24 = strainls;

load('linear-elasticity/tip_variety_data_pois_0.3_young_30.mat');
a_l = a;
b_l = b;
gam_peaks_l = gam_peaks;
R0_l = R0;
RN_l = RN;
tip_lengths_l = tip_lengths;
ZN_l = ZN;
W_ints_l = W_ints;
W_l = W;
gams_l = gams;
L0Splines_l = L0Splines;
strainls_l = strainls;

load('tip_variety_data_hyphoids_k_6_mu_6.mat');
a_h_nl_6 = a;
b_h_nl_6 = b;
gam_peaks_h_nl_6 = gam_peaks;
R0_h_nl_6 = R0;
RN_h_nl_6 = RN;
tip_lengths_h_nl_6 = tip_lengths;
ZN_h_nl_6 = ZN;
W_ints_h_nl_6 = W_ints;
W_h_nl_6 = W;
gams_h_nl_6 = gams;
L0Splines_h_nl_6 = L0Splines;
strainls_h_nl_6 = strainls;

load('tip_variety_data_hyphoids_k_12_mu_12.mat');
a_h_nl_12 = a;
b_h_nl_12 = b;
gam_peaks_h_nl_12 = gam_peaks;
R0_h_nl_12 = R0;
RN_h_nl_12 = RN;
tip_lengths_h_nl_12 = tip_lengths;
ZN_h_nl_12 = ZN;
W_ints_h_nl_12 = W_ints;
W_h_nl_12 = W;
gams_h_nl_12 = gams;
L0Splines_h_nl_12 = L0Splines;
strainls_h_nl_12 = strainls;

load('tip_variety_data_hyphoids_k_24_mu_24.mat');
a_h_nl_24 = a;
b_h_nl_24 = b;
gam_peaks_h_nl_24 = gam_peaks;
R0_h_nl_24 = R0;
RN_h_nl_24 = RN;
tip_lengths_h_nl_24 = tip_lengths;
ZN_h_nl_24 = ZN;
W_ints_h_nl_24 = W_ints;
W_h_nl_24 = W;
gams_h_nl_24 = gams;
L0Splines_h_nl_24 = L0Splines;
strainls_h_nl_24 = strainls;

load('linear-elasticity/tip_variety_data_hyphoids_pois_0.3_young_30.mat');
a_h_l = a;
b_h_l = b;
gam_peaks_h_l = gam_peaks;
R0_h_l = R0;
RN_h_l = RN;
tip_lengths_h_l = tip_lengths;
ZN_h_l = ZN;
W_ints_h_l = W_ints;
W_h_l = W;
gams_h_l = gams;
L0Splines_h_l = L0Splines;
strainls_h_l = strainls;

%% load data from real-world cell profiles
load('../cell-profiles/dumais-geometry-and-secretion.mat');
gam_peaks_root_hair = gam_peaks;
R0_root_hair = R0;
RN_root_hair = RN;
tip_lengths_root_hair = tip_lengths;
ZN_root_hair = ZN;
W_ints_root_hair = W_ints;
W_root_hair = W;
load('../cell-profiles/rojas-geometry-and-secretion.mat');
gam_peaks_pollen_tube = gam_peaks;
R0_pollen_tube = R0;
RN_pollen_tube = RN;
tip_lengths_pollen_tube = tip_lengths;
ZN_pollen_tube = ZN;
W_ints_pollen_tube = W_ints;
W_pollen_tube = W;

load('../cell-profiles/abenza-yeast-canonical.mat');
ZN_yeast = [ZN];
RN_yeast = [RN];
R0_yeast = [Rtip];
tip_lengths_yeast = [tip_length];
gam_peaks_yeast = [gam_peak];
W_int_yeast = W_int;
W_yeast = [W_int(2) - W_int(1)];

% 3D plot - s_p/L
hold on;
plot3(ZN_nl_6 ./ RN_nl_6, R0_nl_6 ./ RN_nl_6, gam_peaks_nl_6 ./ tip_lengths_nl_6, 'LineWidth', 2.0); % secretion peaks
plot3(ZN_nl_12 ./ RN_nl_12, R0_nl_12 ./ RN_nl_12, gam_peaks_nl_12 ./ tip_lengths_nl_12, 'LineWidth', 2.0);
plot3(ZN_nl_24 ./ RN_nl_24, R0_nl_24 ./ RN_nl_24, gam_peaks_nl_24 ./ tip_lengths_nl_24, 'LineWidth', 2.0);
plot3(ZN_l ./ RN_l, R0_l ./ RN_l, gam_peaks_l ./ tip_lengths_l, 'LineWidth', 2.0);
scatter3(ZN_root_hair./RN_root_hair, R0_root_hair./RN_root_hair, gam_peaks_root_hair./tip_lengths_root_hair, 'LineWidth', 2.0);
scatter3(ZN_pollen_tube./RN_pollen_tube, R0_pollen_tube./RN_pollen_tube, gam_peaks_pollen_tube./tip_lengths_pollen_tube, 'LineWidth', 2.0);
scatter3(ZN_yeast./RN_yeast, R0_yeast./RN_yeast, gam_peaks_yeast./tip_lengths_yeast, 'LineWidth', 2.0);
grid on;
xlabel('ZN/RN');
% xlim([0.7 1.1]);
ylabel('R0/RN');
% ylim([0.9 1.3]);
zlabel('s_p/L');
zlim([0 1]);
legend('Nonlinear: K = \mu = 6', 'Nonlinear: K = \mu = 12', 'Nonlinear: K = \mu = 24', 'Linear: \nu = 0.3, E = 30', 'Root hair', 'Pollen tube', 'Yeast (canonical data)');
saveas(gcf, 'tip_shape_variety_gam_peak_3d_aggregate_plot.fig');
close all;

% 2D plot - s_p/L
hold on;
plot(R0_nl_6 ./ RN_nl_6, gam_peaks_nl_6 ./ tip_lengths_nl_6, 'LineWidth', 2.0); % secretion peaks - ellipsoids
plot(R0_nl_12 ./ RN_nl_12, gam_peaks_nl_12 ./ tip_lengths_nl_12, 'LineWidth', 2.0);
plot(R0_nl_24 ./ RN_nl_24, gam_peaks_nl_24 ./ tip_lengths_nl_24, 'LineWidth', 2.0);
plot(R0_l ./ RN_l, gam_peaks_l ./ tip_lengths_l, 'LineWidth', 2.0);
plot(R0_h_nl_6 ./ RN_h_nl_6, gam_peaks_h_nl_6 ./ tip_lengths_h_nl_6, 'LineWidth', 2.0); % secretion peaks - hyphoids
plot(R0_h_nl_12 ./ RN_h_nl_12, gam_peaks_h_nl_12 ./ tip_lengths_h_nl_12, 'LineWidth', 2.0);
plot(R0_h_nl_24 ./ RN_h_nl_24, gam_peaks_h_nl_24 ./ tip_lengths_h_nl_24, 'LineWidth', 2.0);
plot(R0_h_l ./ RN_h_l, gam_peaks_h_l ./ tip_lengths_h_l, 'LineWidth', 2.0);
scatter(R0_root_hair./RN_root_hair, gam_peaks_root_hair./tip_lengths_root_hair, 'LineWidth', 2.0);
scatter(R0_pollen_tube./RN_pollen_tube, gam_peaks_pollen_tube./tip_lengths_pollen_tube, 'LineWidth', 2.0);
grid on;
xlabel('R_0/R_N');
xlim([0.57 1.35]);
ylabel('s_p/L');
ylim([0 1]);
l = legend('Ellipsoid: nonlinear K = \mu = 6', 'Ellipsoid: nonlinear K = \mu = 12', 'Ellipsoid: nonlinear K = \mu = 24', 'Ellipsoid: linear \nu = 0.3, E = 30', ...
    'Hyphoid: nonlinear K = \mu = 6', 'Hyphoid: nonlinear K = \mu = 12', 'Hyphoid: nonlinear K = \mu = 24', 'Hyphoid: linear \nu = 0.3, E = 30', ...
    'Root hair', 'Pollen tube', 'Location', 'northwest');
set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.8]));
set(gca, 'fontsize', 12);
exportgraphics(gcf, 'tip_shape_variety_gam_peak_with_data.png');
close all;

% 3D plot - W/L
hold on;
plot3(ZN_nl_6 ./ RN_nl_6, R0_nl_6 ./ RN_nl_6, W_nl_6 ./ tip_lengths_nl_6, 'LineWidth', 2.0);
plot3(ZN_nl_12 ./ RN_nl_12, R0_nl_12 ./ RN_nl_12, W_nl_12 ./ tip_lengths_nl_12, 'LineWidth', 2.0);
plot3(ZN_nl_24 ./ RN_nl_24, R0_nl_24 ./ RN_nl_24, W_nl_24 ./ tip_lengths_nl_24, 'LineWidth', 2.0);
plot3(ZN_l ./ RN_l, R0_l ./ RN_l, W_l ./ tip_lengths_l, 'LineWidth', 2.0);
scatter3(ZN_root_hair./RN_root_hair, R0_root_hair./RN_root_hair, W_root_hair./tip_lengths_root_hair, 'LineWidth', 2.0);
scatter3(ZN_pollen_tube./RN_pollen_tube, R0_pollen_tube./RN_pollen_tube, W_pollen_tube./tip_lengths_pollen_tube, 'LineWidth', 2.0);
scatter3(ZN_yeast./RN_yeast, R0_yeast./RN_yeast, W_yeast./tip_lengths_yeast, 'LineWidth', 2.0);
grid on;
xlabel('ZN/RN');
% xlim([0.7 1.1]);
ylabel('R0/RN');
% ylim([0.9 1.3]);
zlabel('W/L');
zlim([0 1]);
legend('Nonlinear: K = \mu = 6', 'Nonlinear: K = \mu = 12', 'Nonlinear: K = \mu = 24', 'Linear: \nu = 0.3, E = 30', 'Root hair', 'Pollen tube', 'Yeast (canonical data)');
saveas(gcf, 'tip_shape_variety_secretion_polarity_3d_aggregate_plot.fig');
close all;

% 2D plot - W/L
hold on;
plot(R0_nl_6 ./ RN_nl_6, W_nl_6 ./ tip_lengths_nl_6, 'LineWidth', 2.0); % secretion peaks
plot(R0_nl_12 ./ RN_nl_12, W_nl_12 ./ tip_lengths_nl_12, 'LineWidth', 2.0);
plot(R0_nl_24 ./ RN_nl_24, W_nl_24 ./ tip_lengths_nl_24, 'LineWidth', 2.0);
plot(R0_l ./ RN_l, W_l ./ tip_lengths_l, 'LineWidth', 2.0);
plot(R0_h_nl_6 ./ RN_h_nl_6, W_h_nl_6 ./ tip_lengths_h_nl_6, 'LineWidth', 2.0); % secretion peaks - hyphoids
plot(R0_h_nl_12 ./ RN_h_nl_12, W_h_nl_12 ./ tip_lengths_h_nl_12, 'LineWidth', 2.0);
plot(R0_h_nl_24 ./ RN_h_nl_24, W_h_nl_24 ./ tip_lengths_h_nl_24, 'LineWidth', 2.0);
plot(R0_h_l ./ RN_h_l, W_h_l ./ tip_lengths_h_l, 'LineWidth', 2.0);
scatter(R0_root_hair./RN_root_hair, W_root_hair./tip_lengths_root_hair, 'LineWidth', 2.0);
scatter(R0_pollen_tube./RN_pollen_tube, W_pollen_tube./tip_lengths_pollen_tube, 'LineWidth', 2.0);
grid on;
xlabel('R_0/R_N');
xlim([0.57 1.35]);
ylabel('W/L');
ylim([0 1]);
l = legend('Ellipsoid: nonlinear K = \mu = 6', 'Ellipsoid: nonlinear K = \mu = 12', 'Ellipsoid: nonlinear K = \mu = 24', 'Ellipsoid: linear \nu = 0.3, E = 30', ...
    'Hyphoid: nonlinear K = \mu = 6', 'Hyphoid: nonlinear K = \mu = 12', 'Hyphoid: nonlinear K = \mu = 24', 'Hyphoid: linear \nu = 0.3, E = 30', ...
    'Root hair', 'Pollen tube', 'Location', 'northwest');
set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.8]));
set(gca, 'fontsize', 12);
exportgraphics(gcf, 'tip_shape_variety_secretion_polarity_with_data.png');
close all;