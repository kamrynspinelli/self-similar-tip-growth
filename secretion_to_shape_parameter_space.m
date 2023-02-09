means = 0:0.1:0.7;
stds = 0.25:0.1:0.65;

% for i = 1:size(means, 2)
%     for j = 1:size(stds, 2)
%         [tr, ev, sr] = secretion_to_shape_search_evolution(means(i), stds(j), 48, 600);
%         tip_radius(i,j) = tr;
%         elong_vel(i,j) = ev;
%         steady_radius(i,j) = sr;
%     end
% end
% 
% save('secretion_to_shape_parameter_space_data.mat');
load('secretion_to_shape_parameter_space_data.mat');
close all;

% heatmap of tip radius of curvature (not normalized)
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, tip_radius);
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$R_0$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_tip_radius_only.png');

% heatmap of steady state radius (not normalized)
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, steady_radius);
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$R_N$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_steady_radius_only.png');

% heatmap plot of log tip radius of curvature / steady state radius
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, log(tip_radius ./ steady_radius));
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$\log (R_0 / R_N)$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_tip_radius.png');

% heatmap plot of nonlog tip radius of curvature / steady state radius
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, tip_radius ./ steady_radius);
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$R_0 / R_N$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_tip_radius_nonlog.png');

% heatmap plot of log elongation velocity / steady state radius
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, log(elong_vel ./ steady_radius));
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$\log (v / (\bar{\gamma} R_N))$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_elongation_velocity.png');

% heatmap plot of nonlog elongation velocity / steady state radius
xlabels = {'0.25', '0.35', '0.45', '0.55', '0.65'};
ylabels = {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'};
h = heatmap(xlabels, ylabels, elong_vel ./ steady_radius);
h.XLabel = '\sigma';
h.YLabel = '\mu';
title('$v / (\bar{\gamma} R_N)$');
h.NodeChildren(3).Title.Interpreter = 'latex'; % hack to use the latex interpreter for title
colormap(turbo);
exportgraphics(gcf, 'media/heatmap_elongation_velocity_nonlog.png');