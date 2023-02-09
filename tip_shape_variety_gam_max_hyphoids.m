function [] = tip_shape_variety_gam_max(stiffness)
min_param = 3.4;
max_param = 9;
param_step = 0.2;

a = min_param:param_step:max_param;

gam_peaks = zeros(size(a));
tip_lengths = zeros(size(a));
ZN = zeros(size(a));
RN = zeros(size(a));
R0 = zeros(size(a));
W_ints = zeros(size(a,2), 2);
W = zeros(size(a));
ext_peaks = zeros(size(a));
W_ints_ext = zeros(size(a,2), 2);
W_ext = zeros(size(a));
ks_peaks = zeros(size(a));

for i = 1:size(a,2)
    [tip_length gam_peak Z R Rtip W_int ext_peak W_int_ext gam L0Spline strainl ksspline] = growth_anisotropic_hyphoid(a(i), stiffness);
    gam_peaks(i) = gam_peak;
    tip_lengths(i) = tip_length;
    ZN(i) = Z;
    RN(i) = R;
    R0(i) = Rtip;
    W_ints(i,:) = W_int;
    W(i) = W_int(2) - W_int(1);
    ext_peaks(i) = ext_peak;
    W_ints_ext(i,:) = W_int_ext;
    W_ext(i) = W_int_ext(2) - W_int_ext(1);
    gams(i,:) = gam;
    L0Splines(i,:) = L0Spline;
    strainls(i,:) = strainl;
    kssplines(i,:) = ksspline;
end
save(['tip_variety_data_hyphoids_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat'], 'gam_peaks', 'tip_lengths', 'a', 'ZN', 'RN', 'R0', 'W_ints', 'W', 'ext_peaks', 'W_ext', 'W_ints_ext', 'gams', 'L0Splines', 'strainls', 'kssplines');

hold on;
fill([a fliplr(a)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1]); % periwinkle fill
plot(a, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min_param max_param]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_hyphoids_parameter_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;

hold on;
fill([ZN./RN fliplr(ZN./RN)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1])
plot(ZN ./ RN, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min(ZN ./ RN) max(ZN ./ RN)]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_hyphoids_cartesian_ratio_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;

hold on;
fill([R0./RN fliplr(R0./RN)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1])
plot(R0 ./ RN, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min(R0 ./ RN) max(R0 ./ RN)]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_hyphoids_radius_ratio_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;
end