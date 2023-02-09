clear all;

means = [0];
% means = [0 0.5 1];
stds = [0.25 0.5];

tip_radius6 = zeros(size(means,2), size(stds,2));
tip_radius24 = zeros(size(means,2), size(stds,2));
elong_vel6 = zeros(size(means,2), size(stds,2));
elong_vel24 = zeros(size(means,2), size(stds,2));
gam_mass6 = zeros(size(means,2), size(stds,2));
gam_mass24 = zeros(size(means,2), size(stds,2));
steady_radius6 = zeros(size(means,2), size(stds,2));
steady_radius24 = zeros(size(means,2), size(stds,2));
for i = 1:size(means, 2)
    for j = 1:size(stds,2)
        [tr24, ev24, gm24, sr24, ss24] = secretion_to_shape_search_evolution(means(i) * stds(j), stds(j), 24/stds(j), 500, 24 * stds(j));
        % N=48 patches for \sigma=0.5, N=96 patches for \sigma=0.25
        tip_radius24(i,j) = tr24;
        elong_vel24(i,j) = ev24;
        gam_mass24(i,j) = gm24;
        steady_radius24(i,j) = sr24;
        strain_shank24(i,j) = ss24;
        [tr6, ev6, gm6, sr6, ss6] = secretion_to_shape_search_evolution(means(i) * stds(j), stds(j), 24/stds(j), 250, 6 * stds(j));
        tip_radius6(i,j) = tr6;
        elong_vel6(i,j) = ev6;
        gam_mass6(i,j) = gm6;
        steady_radius6(i,j) = sr6;
        strain_shank6(i,j) = ss6;
    end
end

1;
save('secretion_to_shape_check_length_scale_independence.mat');