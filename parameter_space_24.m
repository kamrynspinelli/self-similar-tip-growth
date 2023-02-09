clear all;
normalized_means = 0:0.1:2.5;

for i = 1:size(normalized_means, 2)
    [tr24, ev24, gm24, sr24, ss24] = secretion_to_shape_search_evolution(normalized_means(i) * 0.25, 0.25, 96, 500, 24 * 0.25);
    tip_radius24(i) = tr24;
    elong_vel24(i) = ev24;
    gam_mass24(i) = gm24;
    steady_radius24(i) = sr24;
    strain_shank24(i) = ss24;
end

save('secretion_to_shape_parameter_space_data_24_truncated.mat');