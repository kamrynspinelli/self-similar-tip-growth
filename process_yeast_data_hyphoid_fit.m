clear all;
load('../cell-profiles/abenza-yeast/actual/cropreg347-01_03_R3D_PRJ-t.mat'); % load the data

N = 64; % the desired number of points for each outline

%% Rotate and process the point data
% compute the rotation angle
spot = 1; % keep track of how many sample points we're using for the least squares computation
for i = 1:size(Cx,1); for j = 1:size(Cx,2) % build the list of points on the tube and old end
    if Cy(i,j) <= 55 / 17 * Cx(i,j) + 127 - 55 * 98 / 17 % && Cy(i,j) >= 55 / 17 * Cx(i,j) + 96 - 55 * 179 / 17
%     if Cy(i,j) <= 55 / 17 * Cx(i,j) + 127 - 75 * 98 / 17
        rot_smple_points(1,spot) = Cx(i,j);
        rot_smple_points(2,spot) = Cy(i,j);
        spot = spot + 1;
    end
end; end

coords = [rot_smple_points(1,:)', rot_smple_points(2,:)']; % SVD magic to find the line minimizing the sum of squared distances: https://stackoverflow.com/a/12251893
[U, S, V] = svd(coords - repmat(mean(coords), size(coords, 1), 1), 0);
a = -V(1, end) / V(2, end);
b = mean(coords * V(:, end)) / V(2, end);

hold on; for i = 1:size(Cx,1); plot(Cx(:,i), Cy(:,i)); end; daspect([1 1 1]);
plot([0 250], a * [0 250] + b);
exportgraphics(gcf, 'media/abenza-yeast-processing-before-rotation.png');
close all;

% do the rotation and save the point cloud
angle = - atan(-1/a) + pi/2;
rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
X = [];
Y = [];
for i = 1:size(Cx,2)
    rottmp = rot * [Cx(:,i) Cy(:,i)-b]';
    Cx_rot(:,i) = rottmp(1,:);
    Cy_rot(:,i) = rottmp(2,:);
    point_indices = find(rottmp(1,:) >= max(rottmp(1,:))-80);
    if 10 <= i && i < 18
        X = [X rottmp(1,point_indices)-max(rottmp(1,point_indices))];
        Y = [Y rottmp(2,point_indices)];
    end
end
hold on; for i = 1:size(Cx,1); plot(Cx_rot(:,i), Cy_rot(:,i)); end; daspect([1 1 1]);
exportgraphics(gcf, 'media/abenza-yeast-processing-after-rotation.png');
close all;

%% Do the fitting
height_scale = (max(Y) - min(Y)) / 2;
X = X / height_scale;
Y = Y / height_scale;

a_hyp = 4:0.1:10;
% a_hyp = 5:0.01:6.5;
for i = 1:size(a_hyp,2)
    score_hyp(i) = 0;
    for j = 1:size(X, 2) % for each point in the point cloud,
%         y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y - Y(j))^2, 0.0001);
%         score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y_closest - Y(j))^2;
        y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (y - Y(j))^2, 0.0001);
        score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (y_closest - Y(j))^2;
    end
end

a_ell = 0.5:0.1:3;
% a_ell = 1.4:0.01:2.4;
for i = 1:size(a_ell,2)
    score_ell(i) = 0;
    for j = 1:size(X, 2) % for each point in the point cloud,
        y_closest = fminunc(@(y) min([(max([0, X(j) - -a_ell(i), -(X(j) - (-3 * a_ell(i)))]))^2 + (Y(j) - 1)^2, ...
            (a_ell(i) * sqrt(1 - y^2) - a_ell(i) - X(j))^2 + (y - Y(j))^2, ...
            (max([0, X(j) - -a_ell(i), -(X(j) - (-8 * a_ell(i)))]))^2 + (Y(j) - -1)^2]), Y(j));
        score_ell(i) = score_ell(i) + min([(max([0, X(j) - -a_ell(i), -(X(j) - (-3 * a_ell(i)))]))^2 + (Y(j) - 1)^2, ...
            (a_ell(i) * sqrt(1 - y_closest^2) - a_ell(i) - X(j))^2 + (y_closest - Y(j))^2, ...
            (max([0, X(j) - -a_ell(i), -(X(j) - (-8 * a_ell(i)))]))^2 + (Y(j) - -1)^2]);
    end
end

fit_hyp = a_hyp(find(score_hyp == min(score_hyp)));
fit_ell = a_ell(find(score_ell == min(score_ell)));

disp(['The best hyphoid approximation for yeast is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(X,2)), ' between each point and the fitted curve.']);
disp(['The best ellipsoid approximation is with a = ', num2str(fit_ell), ', with average squared distance ', num2str(min(score_ell) / size(X,2)), ' between each point and the fitted curve.']);

hold on; plot(X, Y, 'x');
% a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96) / base_height(fit_hyp), 'LineWidth', 2.0); 
a = fit_hyp; plot(pi / a * (-0.963:0.001:0.963) .* cot(pi * (-0.963:0.001:0.963)) - 1/a, (-0.963:0.001:0.963), 'LineWidth', 2.0); 
a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
set(gca, 'ColorOrderIndex', 3);
plot([pi / fit_hyp * 0.963 * cot(pi * 0.963) - 1/fit_hyp, -fit_ell], [1, 1], 'LineWidth', 2.0); % cylinder region for the ellipsoid tip
set(gca, 'ColorOrderIndex', 3);
plot([pi / fit_hyp * 0.963 * cot(pi * 0.963) - 1/fit_hyp, -fit_ell], [-1, -1], 'LineWidth', 2.0);
daspect([1 1 1]);
legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(X, 2))], ...
    ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(X, 2))], 'Location', 'west');
ax = gca;
set(ax, 'fontsize', 12);
xlim tight; ylim([-1.1 1.1]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4]);
pause(0.5);
exportgraphics(gcf, ['media/yeast_hyphoid_ellpisoid_fit.png']);
close all;