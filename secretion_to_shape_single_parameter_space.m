normalized_means = 0:0.1:3;
sigma = 0.25; % length scale used for the forward simulations

% for i = 1:size(normalized_means, 2)
%     [tr24, ev24, sr24] = secretion_to_shape_search_evolution(normalized_means(i) / 2, 0.5, 48, 600, 24);
%     tip_radius24(i) = tr24;
%     elong_vel24(i) = ev24;
%     steady_radius24(i) = sr24;
%     [tr6, ev6, sr6] = secretion_to_shape_search_evolution(normalized_means(i) / 2, 0.5, 48, 600, 6);
%     tip_radius6(i) = tr6;
%     elong_vel6(i) = ev6;
%     steady_radius6(i) = sr6;
% end
% 
% save('secretion_to_shape_parameter_space_data.mat');
% close all;

% load('secretion_to_shape_parameter_space_data.mat');

% load the data from the forward direction
% load('secretion_to_shape_parameter_space_data_6_incomplete.mat');
% load('secretion_to_shape_parameter_space_data_24.mat');
load('secretion_to_shape_parameter_space_data_24_with_shank_r0.mat');

% compute l and W for each of the prescribed distributions
for i = 1:size(normalized_means, 2)
    mean = normalized_means(i);
    gam = @(x) 1/(sqrt(2*pi)) * (exp(-1/2 * (x-mean).^2) + exp(-1/2 * (x+mean).^2));
    
     % find the approximate location of the maximum
    gam_num = gam(0:0.001:(4*max(normalized_means)));
    ell_gaussian(i) = 0.001 * (max(find(gam_num == max(gam_num))) - 1);
    
    levels = fliplr(unique(gam_num));
    s1_ind = find(gam_num == max(gam_num));
    s2_ind = find(gam_num == max(gam_num));
    total_gam = 1;
    window_int = sum(gam_num(s1_ind:s2_ind-1)) * 0.001;
    l = 1;
    while window_int < 0.5 * total_gam % assume that \gamma is either monotonic or has a single peak
        l = l + 1;
        p = find(gam_num == levels(l));
        if p < s1_ind
            s1_ind = p;
        else
            s2_ind = p;
        end
        window_int = sum(gam_num(s1_ind:s2_ind-1)) * 0.001;
    end
    s1 = 0.001 * (s1_ind - 1);
    s2 = 0.001 * (s2_ind - 1);
    W_int = [s1 s2];
    if i == 26
        1;
    end
    W_gaussian(i) = s2-s1;
end
W_gaussian = sigma * W_gaussian; % correction for length scale sigma
ell_gaussian = sigma * ell_gaussian;

% load the data from the inverse direction
load('tip_variety_data_hyphoids_k_6_mu_6.mat');
a_h_nl_6 = a;
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
strain_shanks_h_nl_6 = strainls(:,1)';
kssplines_h_nl_6 = kssplines;
max_ks_h_nl_6 = max(kssplines, [], 2);
tip_ks_h_nl_6 = kssplines(:,end);
max_gams_h_nl_6 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_nl_6(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
    elong_vel_h_nl_6(i) = elongation_velocity(gams(i,:), strainls(i,:)', L0Splines(i,:));
    gam_mass_h_nl_6(i) = gam_mass(gams(i,:), strainls(i,:)', L0Splines(i,:));
end

load('tip_variety_data_hyphoids_k_24_mu_24.mat');
a_h_nl_24 = a;
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
strain_shanks_h_nl_24 = strainls(:,1)';
kssplines_h_nl_24 = kssplines;
max_ks_h_nl_24 = max(kssplines, [], 2);
tip_ks_h_nl_24 = kssplines(:,end);
max_gams_h_nl_24 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_nl_24(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
    elong_vel_h_nl_24(i) = elongation_velocity(gams(i,:), strainls(i,:)', L0Splines(i,:));
    gam_mass_h_nl_24(i) = gam_mass(gams(i,:), strainls(i,:)', L0Splines(i,:));
end

% R0RNmin = min([(tip_radius6 ./ steady_radius6) (tip_radius24 ./ steady_radius24)]);
% R0RNmax = max([(tip_radius6 ./ steady_radius6) (tip_radius24 ./ steady_radius24)]);

hold on;
% plot(normalized_means(1:26), tip_radius6(1:26), normalized_means(1:26), tip_radius24(1:26), 'LineWidth', 2.0);
plot(normalized_means(1:26), tip_radius24(1:26) / sigma, 'LineWidth', 2.0);
xlabel('$\mu / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$R_0 / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0 2.5]);
% legend('$\tilde{K}_h = \tilde{\mu}_h = 6$', '$\tilde{K}_h = \tilde{\mu}_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'media/secretion_to_shape_tip_radius_only_nrmlz_mean.png');
close all;

hold on;
% plot(normalized_means(1:26), steady_radius6(1:26), normalized_means(1:26), steady_radius24(1:26), 'LineWidth', 2.0);
plot(normalized_means(1:26), steady_radius24(1:26) / sigma, 'LineWidth', 2.0);
xlabel('$\mu / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$R_N / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0 2.5]);
% legend('$K_h = \mu_h = 6$', '$K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'media/secretion_to_shape_steady_radius_only_nrmlz_mean.png');
close all;

hold on;
% plot(normalized_means(1:26), tip_radius6(1:26) ./ steady_radius6(1:26), normalized_means(1:26), tip_radius24(1:26) ./ steady_radius24(1:26), 'LineWidth', 2.0);
plot(normalized_means(1:26), tip_radius24(1:26) ./ steady_radius24(1:26), 'LineWidth', 2.0);
xlabel('$\mu / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$R_0 / R_N$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0 2.5]);
% legend('$K_h = \mu_h = 6$', '$K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'media/secretion_to_shape_tip_radius_nrmlz_mean.png');
close all;

% hold on;
% plot(normalized_means(1:26), elong_vel6(1:26), normalized_means(1:26), elong_vel24(1:26), 'LineWidth', 2.0);
% xlabel('$\mu / \sigma$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$v / V$', 'Interpreter', 'latex', 'FontSize', 14);
% xlim([0 3]);
% % legend('$K_h = \mu_h = 6$', '$K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'east');
% exportgraphics(gcf, 'media/secretion_to_shape_elongation_velocity_nrmlz_mean.png');
% close all;

hold on;
% plot(tip_radius6(1:26) ./ steady_radius6(1:26), ell_gaussian(1:26) ./ steady_radius6(1:26), tip_radius24(1:26) ./ steady_radius24(1:26), ell_gaussian(1:26) ./ steady_radius24(1:26), 'LineWidth', 2.0);
plot(tip_radius24(1:26) / sigma, ell_gaussian(1:26) / sigma, 'LineWidth', 2.0);
% plot(R0_h_nl_6 ./ RN_h_nl_6, gam_peaks_h_nl_6 ./ RN_h_nl_6, R0_h_nl_24 ./ RN_h_nl_24, gam_peaks_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 2.0);
xlabel('$\tilde{R}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\tilde{\ell}$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'fontsize', 12);
% xlim([R0RNmin R0RNmax]);
% xlim([0.5 1.4]);
% xticks([0.5 0.8 1.1 1.4]);
% ylim([0 2.5]);
% yticks([0 0.5 1 1.5 2 2.5]);
% xlim([1 10]);
% xticks([1 4 7 10]);
xlim([1 9.5]);
xticks([2 4 6 8]);
ylim([0 2.5]);
yticks([0 0.5 1 1.5 2 2.5]);
% legend('$\tilde{K}_h = \tilde{\mu}_h = 24$', 'reverse $K_h = \mu_h = 6$', 'reverse $K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'southeast');
set(gca, 'FontSize', 16);
% inset plots of example tip shapes
ax = axes('Position',[.26 .17 .1 .1]); % inset plot 1: pointy tip
load('forward_direction_inset_1.mat');
plot(rv(1,148:end), rv(2,148:end), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
ax = axes('Position',[.38 .51 .1 .1]); % inset plot 2: medium tip
load('forward_direction_inset_2.mat');
plot(rv(1,180:end), rv(2,180:end), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
ax = axes('Position',[.77 .77 .1 .1]); % inset plot 2: flat tip
load('forward_direction_inset_3.mat');
plot(rv(1,212:end), rv(2,212:end), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
exportgraphics(gcf, 'media/secretion_to_shape_ell_tip_radius.png');
close all;

hold on;
plot(tip_radius24(1:26) / sigma, W_gaussian(1:26) / sigma, 'LineWidth', 2.0);
% plot(R0_h_nl_6 ./ RN_h_nl_6, W_h_nl_6 ./ RN_h_nl_6, R0_h_nl_24 ./ RN_h_nl_24, W_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 2.0);
xlabel('$\tilde{R}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\tilde{w}$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'fontsize', 12);
% xlim([R0RNmin R0RNmax]);
% xlim([0.5 1.4]);
% xticks([0.5 0.8 1.1 1.4]);
% ylim([0 2.5]);
% yticks([0 0.5 1 1.5 2 2.5]);
% xlim([1 10]);
% xticks([1 4 7 10]);
xlim([1 9.5]);
xticks([2 4 6 8]);
ylim([0.6 1.4]);
yticks([0.6 0.8 1 1.2 1.4]);
% legend('$K_h = \mu_h = 6$', '$K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'media/secretion_to_shape_W_tip_radius.png');
close all;

hold on;
% plot(tip_radius24(1:26) / sigma, elong_vel24(1:26) ./ gam_mass24(1:26) ./ (strain_shank24(1:26)-1), 'LineWidth', 2.0);
plot(tip_radius24(1:26) / sigma, elong_vel24(1:26) ./ gam_mass24(1:26) * 2 * 2 * 6 ./ r0_shank24, 'LineWidth', 2.0);
% plot(R0_h_nl_6 ./ RN_h_nl_6, elong_vel_h_nl_6 ./ gam_mass_h_nl_6 ./ (strain_shanks_h_nl_6-1), R0_h_nl_24 ./ RN_h_nl_24, elong_vel_h_nl_24 ./ gam_mass_h_nl_24 ./ (strain_shanks_h_nl_24-1), 'LineWidth', 2.0);
xlabel('$R$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$2 v K_h / V P$', 'Interpreter', 'latex', 'FontSize', 14);
% xlim([min([(R0_h_nl_6 ./ RN_h_nl_6) (R0_h_nl_24 ./ RN_h_nl_24)]) max([(R0_h_nl_6 ./ RN_h_nl_6) (R0_h_nl_24 ./ RN_h_nl_24)])]);
set(gca, 'fontsize', 12);
% xticks([0.6 0.8 1 1.2]);
% xlim([0.5 1.4]);
% xticks([0.5 0.8 1.1 1.4]);
% ylim([0.3 1.2]);
% yticks([0.3 0.6 0.9 1.2]);
% xlim([1 10]);
% xticks([1 4 7 10]);
xlim([1 9.5]);
xticks([2 4 6 8]);
ylim([0.85 1.3]);
yticks([0.85 1 1.15 1.3]);
% legend('$K_h = \mu_h = 6$', '$K_h = \mu_h = 24$', 'inverse $K_h = \mu_h = 6$', 'inverse $K_h = \mu_h = 24$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'media/secretion_to_shape_elongation_velocity_tip_radius.png');
close all;


function ev = elongation_velocity(gam_num, strainl, L0Spline)
%     velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_num(1:end-1) .* fliplr(strainl(3:end) - 1))]); % velocity relative to the tip; last entry is the tip
%     ev = velocity(1);
    ev = strainl(2) * sum(fliplr(L0Spline) .* gam_num .* fliplr(strainl(2:end)' - 1));
end

function m = gam_mass(gam_num, strainl, L0Spline)
    m = sum(gam_num .* fliplr(L0Spline) .* fliplr(strainl(2:end)'));
end