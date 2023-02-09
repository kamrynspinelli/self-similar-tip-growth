clear all;
normalized_means = 0:0.1:2.5;

for i = 1:size(normalized_means, 2)
%     [tr6, ev6, gm6, sr6, ss6] = secretion_to_shape_search_evolution(normalized_means(i) * 0.25, 0.25, 96, 250, 6 * 0.25);
    [tr6, ev6, gm6, sr6, ss6] = secretion_to_shape_search_evolution(normalized_means(i) * 0.25, 0.25, 96, 500, 6 * 0.25);
    tip_radius6(i) = tr6;
    elong_vel6(i) = ev6;
    gam_mass6(i) = gm6;
    steady_radius6(i) = sr6;
    strain_shank6(i) = ss6;
end

save('secretion_to_shape_parameter_space_data_6_truncated.mat');