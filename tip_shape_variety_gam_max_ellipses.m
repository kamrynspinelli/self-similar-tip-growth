function [] = tip_shape_variety_gam_max(stiffness)
min_ratio = 0.78;
% max_ratio = 0.68;
% min_ratio = 0.7;
% max_ratio = 1.1;
max_ratio = 1.5;
ratio_step = 0.02;

aspects = min_ratio:ratio_step:max_ratio; % tip aspect ratios: a/b
gam_peaks = zeros(size(aspects));
tip_lengths = zeros(size(aspects));
a = zeros(size(aspects));
b = zeros(size(aspects));
ZN = zeros(size(aspects));
RN = zeros(size(aspects));
R0 = zeros(size(aspects));
W_ints = zeros(size(aspects,2), 2);
W = zeros(size(aspects));

for i = 1:size(aspects,2)
    [tip_length gam_peak Z R Rtip W_int gam L0Spline strainl ksspline] = growth_anisotropic_ellipse_spline(aspects(i), 1, stiffness);
    gam_peaks(i) = gam_peak;
    tip_lengths(i) = tip_length;
    a(i) = aspects(i);
    b(i) = 1;
    ZN(i) = Z;
    RN(i) = R;
    R0(i) = Rtip;
    W_ints(i,:) = W_int;
    W(i) = W_int(2) - W_int(1);
    gams(i,:) = gam;
    L0Splines(i,:) = L0Spline;
    strainls(i,:) = strainl;
    kssplines(i,:) = ksspline;
end
save(['tip_variety_data_ellipses_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat'], 'gam_peaks', 'tip_lengths', 'a', 'b', 'ZN', 'RN', 'R0', 'W_ints', 'W', 'gams', 'L0Splines', 'strainls', 'kssplines');

hold on;
fill([aspects fliplr(aspects)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1]); % periwinkle fill
plot(aspects, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min_ratio max_ratio]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_ellipses_aspect_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;

hold on;
fill([ZN./RN fliplr(ZN./RN)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1])
plot(ZN ./ RN, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min(ZN ./ RN) max(ZN ./ RN)]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_ellipses_cartesian_ratio_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;

hold on;
fill([R0./RN fliplr(R0./RN)], [W_ints(:,2)'./tip_lengths fliplr(W_ints(:,1)'./tip_lengths)], [0.5 0.6 1])
plot(R0 ./ RN, gam_peaks ./ tip_lengths, 'LineWidth', 2.0);
xlim([min(R0 ./ RN) max(R0 ./ RN)]);
ylim([0 1]);
pbaspect([1 1 1]);
exportgraphics(gcf, ['media/anisotropic_growth_ellipses_radius_ratio_gam_max_stiff_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.png']);
close all;
end