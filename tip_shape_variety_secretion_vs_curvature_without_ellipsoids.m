%% load the ellipsoid and hyphoid trendlines
load('tip_variety_data_ellipses_k_6_mu_6.mat');
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
kssplines_nl_6 = kssplines;
max_ks_nl_6 = max(kssplines, [], 2);
tip_ks_nl_6 = kssplines(:,end);
max_gams_nl_6 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_nl_6(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

load('tip_variety_data_ellipses_k_12_mu_12.mat');
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
kssplines_nl_12 = kssplines;
max_ks_nl_12 = max(kssplines, [], 2);
tip_ks_nl_12 = kssplines(:,end);
max_gams_nl_12 = max(gams, [], 2);
for i = 1:size(gams, 1)
    % find the location of the peak in curvature
    ks_peaks_nl_12(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
    % TESTING
    % find the arclength window containing half of the total curvature
    ks = kssplines(i,:);
    levels = fliplr(unique(ks));
    s1_ind = find(ks == max(ks));
    s2_ind = find(ks == max(ks));
    total_ks = pi/2; % a simple closed plane curve must have total curvature 2pi
    window_int = sum(ks(s1_ind:s2_ind-1) .* fliplr(L0Splines(i,end-s2_ind+2:end-s1_ind+1)) .* fliplr(strainls(i,end-s2_ind+2:end-s1_ind+1)));
    l = 1;
    while window_int < 0.5 * total_ks % assume that \gamma is either monotonic or has a single peak
        l = l + 1;
        p = find(ks == levels(l));
        if p < s1_ind
            s1_ind = p;
        else
            s2_ind = p;
        end
        window_int = sum(ks(s1_ind:s2_ind-1) .* fliplr(L0Splines(i,end-s2_ind+2:end-s1_ind+1)) .* fliplr(strainls(i,end-s2_ind+2:end-s1_ind+1)));
    end
    s1 = sum(fliplr(L0Splines(i,end-s1_ind+2:end)) .* fliplr(strainls(i,end-s1_ind+2:end)));
    s2 = sum(fliplr(L0Splines(i,end-s2_ind+2:end)) .* fliplr(strainls(i,end-s2_ind+2:end)));
%     KW = [s1 s2];
    KW(i) = s2 - s1;
end

load('tip_variety_data_ellipses_k_24_mu_24.mat');
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
kssplines_nl_24 = kssplines;
max_ks_nl_24 = max(kssplines, [], 2);
tip_ks_nl_24 = kssplines(:,end);
max_gams_nl_24 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_nl_24(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

% load('linear-elasticity/tip_variety_data_pois_0.3_young_30.mat');
% a_l = a;
% b_l = b;
% gam_peaks_l = gam_peaks;
% R0_l = R0;
% RN_l = RN;
% tip_lengths_l = tip_lengths;
% ZN_l = ZN;
% W_ints_l = W_ints;
% W_l = W;
% gams_l = gams;
% L0Splines_l = L0Splines;
% strainls_l = strainls;

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
kssplines_h_nl_6 = kssplines;
max_ks_h_nl_6 = max(kssplines, [], 2);
tip_ks_h_nl_6 = kssplines(:,end);
max_gams_h_nl_6 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_nl_6(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
    elong_vel_h_nl_6(i) = elongation_velocity(gams(i,:), strainls(i,:)', L0Splines(i,:));
    gam_mass_h_nl_6(i) = gam_mass(gams(i,:), strainls(i,:)', L0Splines(i,:));
end

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
kssplines_h_nl_12 = kssplines;
max_ks_h_nl_12 = max(kssplines, [], 2);
tip_ks_h_nl_12 = kssplines(:,end);
max_gams_h_nl_12 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_nl_12(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

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
kssplines_h_nl_24 = kssplines;
max_ks_h_nl_24 = max(kssplines, [], 2);
tip_ks_h_nl_24 = kssplines(:,end);
max_gams_h_nl_24 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_nl_24(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';    
    elong_vel_h_nl_24(i) = elongation_velocity(gams(i,:), strainls(i,:)', L0Splines(i,:));
    gam_mass_h_nl_24(i) = gam_mass(gams(i,:), strainls(i,:)', L0Splines(i,:));
end

load('linear-elasticity/tip_variety_data_hyphoids_pois_0.3_young_30.mat');
a_h_l_03_30 = a;
b_h_l_03_30 = b;
gam_peaks_h_l_03_30 = gam_peaks;
R0_h_l_03_30 = R0;
RN_h_l_03_30 = RN;
tip_lengths_h_l_03_30 = tip_lengths;
ZN_h_l_03_30 = ZN;
W_ints_h_l_03_30 = W_ints;
W_h_l_03_30 = W;
gams_h_l_03_30 = gams;
L0Splines_h_l_03_30 = L0Splines;
strainls_h_l_03_30 = strainls;
kssplines_h_l_03_30 = kssplines;
max_ks_h_l_03_30 = max(kssplines, [], 2);
tip_ks_h_l_03_30 = kssplines(:,end);
max_gams_h_l_03_30 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_l_03_30(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

load('linear-elasticity/tip_variety_data_hyphoids_pois_0_young_30.mat');
a_h_l_0_30 = a;
b_h_l_0_30 = b;
gam_peaks_h_l_0_30 = gam_peaks;
R0_h_l_0_30 = R0;
RN_h_l_0_30 = RN;
tip_lengths_h_l_0_30 = tip_lengths;
ZN_h_l_0_30 = ZN;
W_ints_h_l_0_30 = W_ints;
W_h_l_0_30 = W;
gams_h_l_0_30 = gams;
L0Splines_h_l_0_30 = L0Splines;
strainls_h_l_0_30 = strainls;
kssplines_h_l_0_30 = kssplines;
max_ks_h_l_0_30 = max(kssplines, [], 2);
tip_ks_h_l_0_30 = kssplines(:,end);
max_gams_h_l_0_30 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_l_0_30(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

load('linear-elasticity/tip_variety_data_hyphoids_pois_0.3_young_60.mat');
a_h_l_03_60 = a;
b_h_l_03_60 = b;
gam_peaks_h_l_03_60 = gam_peaks;
R0_h_l_03_60 = R0;
RN_h_l_03_60 = RN;
tip_lengths_h_l_03_60 = tip_lengths;
ZN_h_l_03_60 = ZN;
W_ints_h_l_03_60 = W_ints;
W_h_l_03_60 = W;
gams_h_l_03_60 = gams;
L0Splines_h_l_03_60 = L0Splines;
strainls_h_l_03_60 = strainls;
kssplines_h_l_03_60 = kssplines;
max_ks_h_l_03_60 = max(kssplines, [], 2);
tip_ks_h_l_03_60 = kssplines(:,end);
max_gams_h_l_03_60 = max(gams, [], 2);
for i = 1:size(gams, 1)
    ks_peaks_h_l_03_60(i) = strainls(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end) * L0Splines(i,end-find(fliplr(kssplines(i,:)) == max(fliplr(kssplines(i,:))))+2:end)';
end

%% load data from real-world cell profiles
load('../cell-profiles/dumais-geometry-and-secretion_k_24_mu_24.mat');
gam_peaks_root_hair_nl_24 = gam_peaks;
R0_root_hair_nl_24 = R0;
RN_root_hair_nl_24 = RN;
tip_lengths_root_hair_nl_24 = tip_lengths;
ZN_root_hair_nl_24 = ZN;
W_ints_root_hair_nl_24 = W_ints;
W_root_hair_nl_24 = W;
kssplines_root_hair_nl_24 = kssplines;
for i = 1:size(kssplines, 2)
    ksspline = kssplines{i};
    max_ks_root_hair_nl_24(i) = max(ksspline);
    tip_ks_root_hair_nl_24(i) = ksspline(end);
    max_gams_root_hair_nl_24(i) = max(gams{i});
    ks_peaks_root_hair_nl_24(i) = strainls{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Splines{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);
end

load('../cell-profiles/dumais-canonical-hyphoid-shape-k-24-mu-24.mat')
gam_peaks_root_hair_nl_canonical = gam_peak;
R0_root_hair_nl_canonical = Rtip;
RN_root_hair_nl_canonical = RN;
tip_lengths_root_hair_nl_canonical = tip_length;
ZN_root_hair_nl_canonical = ZN;
W_ints_root_hair_nl_canonical = W_int;
W_root_hair_nl_canonical = W_int(2) - W_int(1);
kssplines_root_hair_nl_canonical = ksspline;
max_ks_root_hair_nl_canonical(i) = max(ksspline);
tip_ks_root_hair_nl_canonical(i) = ksspline(end);
max_gams_root_hair_nl_canonical = max(gam);
ks_peaks_root_hair_nl_canonical = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

load('../cell-profiles/dumais-individual-fits-geometry-and-secretion_k_24_mu_24.mat');
gam_peaks_root_hair_individual_fits_nl_24 = gam_peaks;
R0_root_hair_individual_fits_nl_24 = R0;
RN_root_hair_individual_fits_nl_24 = RN;
tip_lengths_root_hair_individual_fits_nl_24 = tip_lengths;
ZN_root_hair_individual_fits_nl_24 = ZN;
W_ints_root_hair_individual_fits_nl_24 = W_ints;
W_root_hair_individual_fits_nl_24 = W;
kssplines_root_hair_individual_fits_nl_24 = kssplines;
for i = 1:size(kssplines, 2)
    ksspline = kssplines{i};
    max_ks_root_hair_individual_fits_nl_24(i) = max(ksspline);
    tip_ks_root_hair_individual_fits_nl_24(i) = ksspline(end);
    max_gams_root_hair_individual_fits_nl_24(i) = max(gams{i});
    ks_peaks_root_hair_individual_fits_nl_24(i) = strainls{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Splines{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);
end

load('../cell-profiles/dumais-geometry-and-secretion_pois_0.3_young_30.mat');
gam_peaks_root_hair_l_03_30 = gam_peaks;
R0_root_hair_l_03_30 = R0;
RN_root_hair_l_03_30 = RN;
tip_lengths_root_hair_l_03_30 = tip_lengths;
ZN_root_hair_l_03_30 = ZN;
W_ints_root_hair_l_03_30 = W_ints;
W_root_hair_l_03_30 = W;
kssplines_root_hair_l_03_30 = kssplines;
for i = 1:size(kssplines, 2)
    ksspline = kssplines{i};
    max_ks_root_hair_l_03_30(i) = max(ksspline);
    tip_ks_root_hair_l_03_30(i) = ksspline(end);
    max_gams_root_hair_l_03_30(i) = max(gams{i});
    ks_peaks_root_hair_l_03_30(i) = strainls{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Splines{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);
end

load('../cell-profiles/dumais-canonical-hyphoid-shape-linelastic-pois-0.3-young-30.mat')
gam_peaks_root_hair_l_canonical = gam_peak;
R0_root_hair_l_canonical = Rtip;
RN_root_hair_l_canonical = RN;
tip_lengths_root_hair_l_canonical = tip_length;
ZN_root_hair_l_canonical = ZN;
W_ints_root_hair_l_canonical = W_int;
W_root_hair_l_canonical = W_int(2) - W_int(1);
kssplines_root_hair_l_canonical = ksspline;
max_ks_root_hair_l_canonical(i) = max(ksspline);
tip_ks_root_hair_l_canonical(i) = ksspline(end);
max_gams_root_hair_l_canonical = max(gam);
ks_peaks_root_hair_l_canonical = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

load('../cell-profiles/dumais-individual-fits-geometry-and-secretion-linelastic_pois_0.3_young_30.mat');
gam_peaks_root_hair_individual_fits_l_03_30 = gam_peaks;
R0_root_hair_individual_fits_l_03_30 = R0;
RN_root_hair_individual_fits_l_03_30 = RN;
tip_lengths_root_hair_individual_fits_l_03_30 = tip_lengths;
ZN_root_hair_individual_fits_l_03_30 = ZN;
W_ints_root_hair_individual_fits_l_03_30 = W_ints;
W_root_hair_individual_fits_l_03_30 = W;
kssplines_root_hair_individual_fits_l_03_30 = kssplines;
for i = 1:size(kssplines, 2)
    ksspline = kssplines{i};
    max_ks_root_hair_individual_fits_l_03_30(i) = max(ksspline);
    tip_ks_root_hair_individual_fits_l_03_30(i) = ksspline(end);
    max_gams_root_hair_individual_fits_l_03_30(i) = max(gams{i});
    ks_peaks_root_hair_individual_fits_l_03_30(i) = strainls{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Splines{i}(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);
end

load('../cell-profiles/yeast-canonical-ellipsoid-shape-k-24-mu-24.mat');
gam_peaks_yeast_nl_canonical = gam_peak;
R0_yeast_nl_canonical = Rtip;
RN_yeast_nl_canonical = RN;
tip_lengths_yeast_nl_canonical = tip_length;
ZN_yeast_nl_canonical = ZN;
W_ints_yeast_nl_canonical = W_int;
W_yeast_nl_canonical = W_int(2) - W_int(1);
kssplines_yeast_nl_canonical = ksspline;
max_ks_yeast_nl_canonical(i) = max(ksspline);
tip_ks_yeast_nl_canonical(i) = ksspline(end);
max_gams_yeast_nl_canonical = max(gam);
ks_peaks_yeast_nl_canonical = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

load('../cell-profiles/moss-caulonema-canonical-hyphoid-shape-k-24-mu-24.mat')
gam_peaks_caulonema_nl_canonical = gam_peak;
R0_caulonema_nl_canonical = Rtip;
RN_caulonema_nl_canonical = RN;
tip_lengths_caulonema_nl_canonical = tip_length;
ZN_caulonema_nl_canonical = ZN;
W_ints_caulonema_nl_canonical = W_int;
W_caulonema_nl_canonical = W_int(2) - W_int(1);
kssplines_caulonema_nl_canonical = ksspline;
max_ks_caulonema_nl_canonical(i) = max(ksspline);
tip_ks_caulonema_nl_canonical(i) = ksspline(end);
max_gams_caulonema_nl_canonical = max(gam);
ks_peaks_caulonema_nl_canonical = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

load('../cell-profiles/moss-chloronema-canonical-hyphoid-shape-k-24-mu-24.mat')
gam_peaks_chloronema_nl_canonical = gam_peak;
R0_chloronema_nl_canonical = Rtip;
RN_chloronema_nl_canonical = RN;
tip_lengths_chloronema_nl_canonical = tip_length;
ZN_chloronema_nl_canonical = ZN;
W_ints_chloronema_nl_canonical = W_int;
W_chloronema_nl_canonical = W_int(2) - W_int(1);
kssplines_chloronema_nl_canonical = ksspline;
max_ks_chloronema_nl_canonical(i) = max(ksspline);
tip_ks_chloronema_nl_canonical(i) = ksspline(end);
max_gams_chloronema_nl_canonical = max(gam);
ks_peaks_chloronema_nl_canonical = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

% load('../cell-profiles/rojas-geometry-and-secretion.mat');
% gam_peaks_pollen_tube = gam_peaks;
% R0_pollen_tube = R0;
% RN_pollen_tube = RN;
% tip_lengths_pollen_tube = tip_lengths;
% ZN_pollen_tube = ZN;
% W_ints_pollen_tube = W_ints;
% W_pollen_tube = W;

% load('../cell-profiles/abenza-yeast-canonical.mat');
% ZN_yeast = [ZN];
% RN_yeast = [RN];
% R0_yeast = [Rtip];
% tip_lengths_yeast = [tip_length];
% gam_peaks_yeast = [gam_peak];
% W_int_yeast = W_int;
% W_yeast = [W_int(2) - W_int(1)];

% % 2D plot - location of max curvature vs location of max secretion (linear
% % only)
% hold on;
% p = plot(ks_peaks_h_l_03_30 ./ RN_h_l_03_30, gam_peaks_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_0_30 ./ RN_h_l_0_30, gam_peaks_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, gam_peaks_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, gam_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
% p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, gam_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% p.MarkerSize = 24; p.Color = 'k';
% q = quiver(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, gam_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, ...
%     (ks_peaks_root_hair_individual_fits_l_03_30 ./ RN_root_hair_individual_fits_l_03_30) - (ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30), ...
%     (gam_peaks_root_hair_individual_fits_l_03_30 ./ RN_root_hair_individual_fits_l_03_30) - (gam_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30), ...
%     0);
% q.MaxHeadSize = 0.15;
% grid on;
% xlabel('s_k / R_N');
% % xlim([0.57 1.35]);
% ylabel('s_p / R_N');
% % ylim([0 1]);
% l = legend('Linear \nu = 0.3, E = 30','Linear \nu = 0, E = 30','Linear \nu = 0.3, E = 60', ...
%     'Root hair (linear \nu = 0.3, E = 30)', 'Root hair (fitted hyphoid, linear)', ...
%     'Location', 'northwest');
% set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% l.FontSize = 8;
% set(gca, 'fontsize', 12);
% daspect([1 1 1]);
% exportgraphics(gcf, 'tip_shape_variety_sk_sp_linear_only.png', 'Resolution', 300);
% close all;

% % 2D plot - location of max curvature vs location of max secretion
% hold on;
% p = plot(ks_peaks_h_nl_6 ./ RN_h_nl_6, gam_peaks_h_nl_6 ./ RN_h_nl_6, 'LineWidth', 1.2); % secretion peaks - hyphoids
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_nl_12 ./ RN_h_nl_12, gam_peaks_h_nl_12 ./ RN_h_nl_12, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_nl_24 ./ RN_h_nl_24, gam_peaks_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_30 ./ RN_h_l_03_30, gam_peaks_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_0_30 ./ RN_h_l_0_30, gam_peaks_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_nl_24 ./ RN_nl_24, gam_peaks_nl_24 ./ RN_nl_24, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, gam_peaks_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% % scatter(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, 'LineWidth', 2.0);
% % scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, gam_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
% p = plot(ks_peaks_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, gam_peaks_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'c');
% % p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, gam_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% % p.MarkerSize = 20; p.Color = 'm';
% % p = plot(ks_peaks_yeast_nl_canonical ./ RN_yeast_nl_canonical, gam_peaks_yeast_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
% p = plot(ks_peaks_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, gam_peaks_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'r');
% % p = plot(ks_peaks_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, gam_peaks_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'g');
% % q = quiver(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, ...
% %     (ks_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     (gam_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     0);
% % q.MaxHeadSize = 0.15;
% grid on;
% xlabel('$\ell_k$', 'Interpreter', 'latex');
% xlim tight;
% xticks(0:0.1:0.7);
% ylabel('$\ell_p$', 'Interpreter', 'latex');
% % ylim tight;
% ylim([0 0.67]);
% l = legend('Nonlinear K_h = \mu_h = 6', 'Nonlinear K_h = \mu_h = 24', ...
%     'Linear \nu = 0.3, E_h = 30','Linear \nu = 0, E_h = 30', ...
%     'Root hair (fitted hyphoid)', 'Multiple species', ...
%     'Location', 'northwest');
% %     'Nonlinear K = \mu = 12', 'Linear \nu = 0.3, E = 60'
% %     'Root hair (nonlinear K = \mu = 24)', 'Root hair (linear \nu = 0.3, E = 30)', ...
% set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% l.FontSize = 7;
% set(gca, 'fontsize', 12);
% pbaspect([1 1 1]);
% h = plot([0.3 0.3], [0 1], '--', 'Color', [0 0 0 0.7], 'LineWidth', 1.5); % vertical divider at s_k/R_n = 0.3 to divide plot
% set(get(get(h,'Annotation'),'LegendInformation'), ... % hack to not display legend entry for the divider
%     'IconDisplayStyle','off');
% % inset plots of example tip shapes
% ax = axes('Position',[.24 .17 .1 .1]); % inset plot 1: pointy tip
% % box off;
% % set(ax, 'Visible','off')
% % set(ax,'YTickLabel',[]);
% a = a_h_nl_24(1);
% h = base_height(a);
% plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
% xlim tight;
% daspect([1 1 1]);
% set(ax,'XTick',[], 'YTick', []);
% ax = axes('Position',[.48 .51 .1 .1]); % inset plot 2: medium tip
% a = a_h_nl_24(9);
% h = base_height(a);
% plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
% xlim tight;
% daspect([1 1 1]);
% set(ax,'XTick',[], 'YTick', []);
% ax = axes('Position',[.71 .6 .1 .1]); % inset plot 3: flat tip
% a = a_h_nl_24(end);
% h = base_height(a);
% plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
% xlim tight;
% daspect([1 1 1]);
% set(ax,'XTick',[], 'YTick', []);
% exportgraphics(gcf, 'tip_shape_variety_sk_sp.png', 'Resolution', 300);
% close all;

% 2D plot - tip radius of curvature vs location of max secretion
hold on;
p = plot(R0_h_nl_6 ./ RN_h_nl_6, gam_peaks_h_nl_6 ./ RN_h_nl_6, 'LineWidth', 1.2); % secretion peaks - hyphoids
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_nl_12 ./ RN_h_nl_12, gam_peaks_h_nl_12 ./ RN_h_nl_12, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_nl_24 ./ RN_h_nl_24, gam_peaks_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(R0_nl_6 ./ RN_nl_6, gam_peaks_nl_6 ./ RN_nl_6, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(R0_nl_24 ./ RN_nl_24, gam_peaks_nl_24 ./ RN_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_l_03_30 ./ RN_h_l_03_30, gam_peaks_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_l_0_30 ./ RN_h_l_0_30, gam_peaks_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(R0_l_0_30 ./ RN_l_0_30, gam_peaks_l_0_30 ./ RN_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_nl_24 ./ RN_nl_24, gam_peaks_nl_24 ./ RN_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, gam_peaks_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% scatter(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, 'LineWidth', 2.0);
% scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, gam_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
p = plot(R0_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, gam_peaks_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'c');
% p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, gam_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% p.MarkerSize = 20; p.Color = 'm';
% p = plot(R0_yeast_nl_canonical ./ RN_yeast_nl_canonical, gam_peaks_yeast_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
p = plot(R0_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, gam_peaks_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
p = plot(R0_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, gam_peaks_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'g');
% q = quiver(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, ...
%     (ks_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
%     (gam_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (gam_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
%     0);
% q.MaxHeadSize = 0.15;
grid on;
xlabel('$\bar{R_0}$', 'Interpreter', 'latex');
xlim tight;
% xticks(0:0.1:0.7);
ylabel('$\ell$', 'Interpreter', 'latex');
% ylim tight;
% ylim([0 0.67]);
l = legend('Nonlinear K_h = \mu_h = 6', 'Nonlinear K_h = \mu_h = 24', ...
    'Linear \nu = 0.3, E_h = 30', 'Linear \nu = 0, E_h = 30', ...
    'Root hair (fitted)', 'Moss chloronema (fitted)', 'Moss caulonema (fitted)', ...
    'Location', 'southeast');
%     'Nonlinear K = \mu = 12', 'Linear \nu = 0.3, E = 60'
%     'Root hair (nonlinear K = \mu = 24)', 'Root hair (linear \nu = 0.3, E = 30)', ...
set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
l.FontSize = 7;
set(gca, 'fontsize', 12);
pbaspect([1 1 1]);
ax = gca;
xsize = ax.XLim(2) - ax.XLim(1);
ysize = ax.YLim(2) - ax.YLim(1);
da = ax.DataAspectRatio;
% h = plot([0.7 0.7], [0 1], '--', 'Color', [0 0 0 0.7], 'LineWidth', 1.5); % vertical divider at s_k/R_n = 0.3 to divide plot
% set(get(get(h,'Annotation'),'LegendInformation'), ... % hack to not display legend entry for the divider
%     'IconDisplayStyle','off');
% inset plots of example tip shapes
ax = axes('Position',[.22 .14 .1 .1]); % inset plot 1: pointy tip
a = a_h_nl_24(1);
h = base_height(a);
plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
ax = axes('Position',[.32 .49 .1 .1]); % inset plot 2: medium tip
a = a_h_nl_24(9);
h = base_height(a);
plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
ax = axes('Position',[.72 .7 .1 .1]); % inset plot 3: flat tip
a = a_h_nl_24(end);
h = base_height(a);
plot(pi / a * (-h:0.001:h) .* cot(pi * (-h:0.001:h)) - 1/a, (-h:0.001:h), 'LineWidth', 2.0)
xlim tight;
daspect([1 1 1]);
set(ax,'XTick',[], 'YTick', []);
exportgraphics(gcf, 'tip_shape_variety_R0bar_lp_no_ellipsoids.png', 'Resolution', 300);
close all;

plot(R0_h_nl_6 ./ RN_h_nl_6, elong_vel_h_nl_6 ./ gam_mass_h_nl_6 * 4 * 6, ...
    R0_h_nl_24 ./ RN_h_nl_24, elong_vel_h_nl_24 ./ gam_mass_h_nl_24 * 4 * 24, 'LineWidth', 8.0)
xlim([0.5 1.25]);
xticks([0.5 1.25]);
ylim([0.55 1.05]);
yticks([0.55 1.05]);
x = xlabel('$R$', 'Interpreter', 'latex');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.03, 0]);
y = ylabel('$V$', 'Interpreter', 'latex');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
ax = gca;
ax.LineWidth = 6.0;
set(gca, 'fontsize', 50);
exportgraphics(gcf, 'tip_shape_variety_elongation_velocity.png', 'Resolution', 300);
close all;

% % 2D plot - location of max curvature vs width of secretion window (linear
% % only)
% hold on;
% p = plot(ks_peaks_h_l_03_30 ./ RN_h_l_03_30, W_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_0_30 ./ RN_h_l_0_30, W_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, W_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, W_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
% p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, W_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% p.MarkerSize = 24; p.Color = 'k';
% q = quiver(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, W_root_hair_l_03_30 ./ RN_root_hair_l_03_30, ...
%     (ks_peaks_root_hair_individual_fits_l_03_30 ./ RN_root_hair_individual_fits_l_03_30) - (ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30), ...
%     (W_root_hair_individual_fits_l_03_30 ./ RN_root_hair_individual_fits_l_03_30) - (W_root_hair_l_03_30 ./ RN_root_hair_l_03_30), ...
%     0);
% q.MaxHeadSize = 0.15;
% grid on;
% xlabel('s_k/R_N');
% % xlim([0 1.4]);
% ylabel('W/R_N');
% % ylim([0 1]);
% l = legend('Linear \nu = 0.3, E = 30','Linear \nu = 0, E = 30','Linear \nu = 0.3, E = 60', ...
%     'Root hair (linear \nu = 0.3, E = 30)', 'Root hair (fitted hyphoid, linear)', ...
%     'Location', 'southwest');
% set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% l.FontSize = 8;
% set(gca, 'fontsize', 12);
% % daspect([1 1 1]);
% exportgraphics(gcf, 'tip_shape_variety_sk_W_linear_only.png', 'Resolution', 300);
% close all;

% % 2D plot - location of max curvature vs width of secretion window
% hold on;
% p = plot(ks_peaks_h_nl_6 ./ RN_h_nl_6, W_h_nl_6 ./ RN_h_nl_6, 'LineWidth', 1.2); % secretion peaks - hyphoids
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_nl_12 ./ RN_h_nl_12, W_h_nl_12 ./ RN_h_nl_12, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_nl_24 ./ RN_h_nl_24, W_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_30 ./ RN_h_l_03_30, W_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_0_30 ./ RN_h_l_0_30, W_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, W_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% % scatter(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, 'LineWidth', 2.0);
% % scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, W_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
% p = plot(ks_peaks_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, W_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'c');
% % p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, W_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% % p.MarkerSize = 20; p.Color = 'm';
% % p = plot(ks_peaks_yeast_nl_canonical ./ RN_yeast_nl_canonical, W_yeast_nl_canonical ./ RN_yeast_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
% % p = plot(ks_peaks_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, W_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'g');
% % p = plot(ks_peaks_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, W_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'r');
% % q = quiver(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, ...
% %     (ks_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     (W_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (W_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     0);
% % q.MaxHeadSize = 0.15;
% grid on;
% xlabel('$\ell_k$', 'Interpreter', 'latex');
% xlim tight;
% xticks(0:0.1:0.7);
% ylabel('$W$', 'Interpreter', 'latex');
% ylim([0.3 0.5]);
% % ylim([0 1]);
% % l = legend('Nonlinear K = \mu = 6', 'Nonlinear K = \mu = 24', ...
% %     'Linear \nu = 0.3, E = 30','Linear \nu = 0, E = 30', ...
% %     'Root hair (fitted hyphoid, nonlinear)', 'Root hair (fitted hyphoid, linear)', ...
% %     'Location', 'south');
% %     'Nonlinear K = \mu = 12', 'Linear \nu = 0.3, E = 60'
% %     'Root hair (nonlinear K = \mu = 24)', 'Root hair (linear \nu = 0.3, E = 30)', ...
% % set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% % l.FontSize = 8;
% set(gca, 'fontsize', 15);
% daspect([1 1 1]);
% h = plot([0.3 0.3], [0 1], '--', 'Color', [0 0 0 0.7], 'LineWidth', 2.0); % vertical divider at s_k/R_n = 0.3 to divide plot
% set(get(get(h,'Annotation'),'LegendInformation'), ... % hack to not display legend entry for the divider
%     'IconDisplayStyle','off');
% exportgraphics(gcf, 'tip_shape_variety_sk_W.png', 'Resolution', 300);
% close all;

% % 2D plot - location of maximum secretion vs width of secretion window
% hold on;
% p = plot(gam_peaks_h_nl_6 ./ RN_h_nl_6, W_h_nl_6 ./ RN_h_nl_6, 'LineWidth', 1.2); % secretion peaks - hyphoids
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_nl_12 ./ RN_h_nl_12, W_h_nl_12 ./ RN_h_nl_12, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% p = plot(gam_peaks_h_nl_24 ./ RN_h_nl_24, W_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(gam_peaks_nl_6 ./ RN_nl_6, W_nl_6 ./ RN_nl_6, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(gam_peaks_nl_24 ./ RN_nl_24, W_nl_24 ./ RN_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(gam_peaks_h_l_03_30 ./ RN_h_l_03_30, W_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(gam_peaks_h_l_0_30 ./ RN_h_l_0_30, W_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% % p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, W_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% % p.Marker = '.'; p.MarkerSize = 10;
% % scatter(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, 'LineWidth', 2.0);
% % scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, W_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
% p = plot(gam_peaks_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, W_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'c');
% % p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, W_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% % p.MarkerSize = 20; p.Color = 'm';
% p = plot(gam_peaks_yeast_nl_canonical ./ RN_yeast_nl_canonical, W_yeast_nl_canonical ./ RN_yeast_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
% p = plot(gam_peaks_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, W_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'r');
% p = plot(gam_peaks_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, W_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'g');
% % q = quiver(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, ...
% %     (ks_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     (W_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (W_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
% %     0);
% % q.MaxHeadSize = 0.15;
% grid on;
% xlabel('$\ell_p$', 'Interpreter', 'latex');
% xlim tight;
% % xticks(0:0.1:0.7);
% ylabel('$W$', 'Interpreter', 'latex');
% % ylim([0.3 0.5]);
% ylim tight;
% % l = legend('Nonlinear K = \mu = 6', 'Nonlinear K = \mu = 24', ...
% %     'Linear \nu = 0.3, E = 30','Linear \nu = 0, E = 30', ...
% %     'Root hair (fitted hyphoid, nonlinear)', 'Root hair (fitted hyphoid, linear)', ...
% %     'Location', 'south');
% %     'Nonlinear K = \mu = 12', 'Linear \nu = 0.3, E = 60'
% %     'Root hair (nonlinear K = \mu = 24)', 'Root hair (linear \nu = 0.3, E = 30)', ...
% % set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% % l.FontSize = 8;
% set(gca, 'fontsize', 15);
% daspect([1 1 1]);
% % h = plot([0.3 0.3], [0 1], '--', 'Color', [0 0 0 0.7], 'LineWidth', 2.0); % vertical divider at s_k/R_n = 0.3 to divide plot
% % set(get(get(h,'Annotation'),'LegendInformation'), ... % hack to not display legend entry for the divider
% %     'IconDisplayStyle','off');
% exportgraphics(gcf, 'tip_shape_variety_lp_W.png', 'Resolution', 300);
% close all;

% 2D plot - tip radius of curvature vs width of secretion window
hold on;
p = plot(R0_h_nl_6 ./ RN_h_nl_6, W_h_nl_6 ./ RN_h_nl_6, 'LineWidth', 1.2); % secretion peaks - hyphoids
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_nl_12 ./ RN_h_nl_12, W_h_nl_12 ./ RN_h_nl_12, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_nl_24 ./ RN_h_nl_24, W_h_nl_24 ./ RN_h_nl_24, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(R0_nl_6 ./ RN_nl_6, W_nl_6 ./ RN_nl_6, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% p = plot(R0_nl_24 ./ RN_nl_24, W_nl_24 ./ RN_nl_24, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_l_03_30 ./ RN_h_l_03_30, W_h_l_03_30 ./ RN_h_l_03_30, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
p = plot(R0_h_l_0_30 ./ RN_h_l_0_30, W_h_l_0_30 ./ RN_h_l_0_30, 'LineWidth', 1.2);
p.Marker = '.'; p.MarkerSize = 10;
% p = plot(ks_peaks_h_l_03_60 ./ RN_h_l_03_60, W_h_l_03_60 ./ RN_h_l_03_60, 'LineWidth', 1.2);
% p.Marker = '.'; p.MarkerSize = 10;
% scatter(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, 'LineWidth', 2.0);
% scatter(ks_peaks_root_hair_l_03_30 ./ RN_root_hair_l_03_30, W_root_hair_l_03_30 ./ RN_root_hair_l_03_30, 'LineWidth', 2.0);
p = plot(R0_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, W_root_hair_nl_canonical ./ RN_root_hair_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'c');
% p = plot(ks_peaks_root_hair_l_canonical ./ RN_root_hair_l_canonical, W_root_hair_l_canonical ./ RN_root_hair_l_canonical, '.', 'LineWidth', 2.0);
% p.MarkerSize = 20; p.Color = 'm';
% p = plot(R0_yeast_nl_canonical ./ RN_yeast_nl_canonical, W_yeast_nl_canonical ./ RN_yeast_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
p = plot(R0_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, W_chloronema_nl_canonical ./ RN_chloronema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'm');
p = plot(R0_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, W_caulonema_nl_canonical ./ RN_caulonema_nl_canonical, '.', 'LineWidth', 2.0, 'MarkerSize', 20, 'Color', 'g');
% q = quiver(ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24, W_root_hair_nl_24 ./ RN_root_hair_nl_24, ...
%     (ks_peaks_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (ks_peaks_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
%     (W_root_hair_individual_fits_nl_24 ./ RN_root_hair_individual_fits_nl_24) - (W_root_hair_nl_24 ./ RN_root_hair_nl_24), ...
%     0);
% q.MaxHeadSize = 0.15;
grid on;
xlabel('$\bar{R}_0$', 'Interpreter', 'latex');
xlim tight;
% xticks(0:0.1:0.7);
ylabel('w', 'Interpreter', 'latex');
% ylim([0.3 0.5]);
ylim tight;
% l = legend('Nonlinear K = \mu = 6', 'Nonlinear K = \mu = 24', ...
%     'Linear \nu = 0.3, E = 30','Linear \nu = 0, E = 30', ...
%     'Root hair (fitted hyphoid, nonlinear)', 'Root hair (fitted hyphoid, linear)', ...
%     'Location', 'south');
%     'Nonlinear K = \mu = 12', 'Linear \nu = 0.3, E = 60'
%     'Root hair (nonlinear K = \mu = 24)', 'Root hair (linear \nu = 0.3, E = 30)', ...
% set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.5]));
% l.FontSize = 8;
ax = gca;
set(gca, 'fontsize', 15);
% daspect([1 1 1]);
daspect(da);
% h = plot([0.3 0.3], [0 1], '--', 'Color', [0 0 0 0.7], 'LineWidth', 2.0); % vertical divider at s_k/R_n = 0.3 to divide plot
% set(get(get(h,'Annotation'),'LegendInformation'), ... % hack to not display legend entry for the divider
%     'IconDisplayStyle','off');
exportgraphics(gcf, 'tip_shape_variety_R0bar_W_no_ellipsoids.png', 'Resolution', 300);
close all;

%% Helper functions
function mid = base_height(a)
    % returns the y s.t. pi*y/a * cot*pi*y) = -2. this is the base of the
    % cell
    f = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
    yl = 0; % bisection method
    yr = 0.5;
    while f(yr) > -2
        yr = (1 + yr) / 2;
    end
    mid = (yl+yr)/2;
    fmid = f(mid);
    while abs(fmid - -2) > 0.001
        if fmid < -2
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f(mid);
    end
end

function ev = elongation_velocity(gam_num, strainl, L0Spline)
%     velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_num(1:end-1) .* fliplr(strainl(3:end) - 1))]); % velocity relative to the tip; last entry is the tip
%     ev = velocity(1);
    ev = strainl(2) * sum(fliplr(L0Spline) .* gam_num .* fliplr(strainl(2:end)' - 1));
end

function m = gam_mass(gam_num, strainl, L0Spline)
    m = sum(gam_num .* fliplr(L0Spline) .* fliplr(strainl(2:end)'));
end