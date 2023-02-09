load('secretion_to_shape_parameter_space_data_24.mat');
normalized_means = 0:0.1:2.5;

for me = 1:size(normalized_means, 2)
    load(['secretion_to_shape_search_evolution_N_96_frames_500_mean_', num2str(0.25*normalized_means(me)), '_std_0.25_stiffness_6.mat']);
    ind = max(find(diff(abs(diff(rv(2,:)) ./ diff(rv(1,:)))) < 0));
    if isempty(ind) % if the first node is the shank
        ind = 1;
    end
    r0_shank24(me) = R0SplineLocal(ind);
end

save('secretion_to_shape_parameter_space_data_24_with_shank_r0.mat');