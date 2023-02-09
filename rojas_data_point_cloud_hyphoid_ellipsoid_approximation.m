clear all;
load('../cell-profiles/PollenOutlines.mat');
indices = [1 5 11 17 21 27 32 37 41 46; ...
    1 5 11 17 21 27 32 37 41 46; ...
    1 5 11 17 21 27 32 37 41 43; ...
    1 5 11 17 21 27 32 37 41 43];
% indices = 1:5:46;
% tipInds = [876 825 820 832 815 757 809 774 780 745]; % indices of the tip points for each outline, computed by finding a point of symmetry for curvature

outline_num = 4;

kssum = zeros(1,60);
kthetasum = zeros(1,60);

Xagg = [];
Yagg = [];

% hold on;
for ind = 1:size(indices, 2)
    X = Outlines_X_Lily{outline_num}(:,indices(outline_num,ind));
    X = X(1:max(find(X)));
    Y = Outlines_Y_Lily{outline_num}(:,indices(outline_num,ind));
    Y = Y(1:max(find(Y)));
    X = X - min(X);
    Y = Y - min(Y);
    
    % find where the lines formed by tangents to the shank intersect
    x1 = X(end);
    y1 = Y(end);
    x2 = X(1);
    y2 = Y(1);
    m1 = (Y(end-1) - Y(end)) / (X(end-1) - X(end));
    m2 = (Y(2) - Y(1)) / (X(2) - X(1));
    intersect =  [1, -m1; 1, -m2] \ [y1 - m1 * x1; y2 - m2 * x2];
    yint = intersect(1);
    xint = intersect(2);
    dists = sqrt((X - xint).^2 + (Y - yint).^2);
    tip_ind = find(dists == min(dists)); % find the point minimizing the distance from the projective point to the shape curve
    m = (yint - Y(tip_ind)) / (xint - X(tip_ind));
    a_hyp = -1/m;
    c = yint/m - xint;
    
    angle = - atan(-1/a_hyp);
    Xtmp = X - X(tip_ind);
    Ytmp = Y - Y(tip_ind);
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    rottmp = rot * [Xtmp Ytmp]';
    Xnew = rottmp(1,:)';
    Ynew = rottmp(2,:)';
%     plot(Xnew, Ynew);
    
    tip_ind = find(Xnew == max(Xnew));
    Xnew = Xnew - Xnew(tip_ind);
    Ynew = Ynew - Ynew(tip_ind);
    Xnew = Xnew / max(Ynew);
    Ynew = Ynew / max(Ynew);
    
    Xagg = [Xagg Xnew(1:40:end)']; % only take some of the points, to avoid waiting hours for the shape optimization
    Yagg = [Yagg Ynew(1:40:end)'];
end

% plot(Xagg, Yagg, 'x'); daspect([1 1 1]);

a_hyp = 2:0.25:12;
% a_hyp = 6.6:0.01:7.6;
for i = 1:size(a_hyp,2)
    score_hyp(i) = 0;
    for j = 1:size(Xagg, 2) % for each point in the point cloud,
        y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - Xagg(j))^2 + (1 / base_height(a_hyp(i)) * y - Yagg(j))^2, 0.0001);
        score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - Xagg(j))^2 + (1 / base_height(a_hyp(i)) * y_closest - Yagg(j))^2;
    end
end

a_ell = 0.5:0.1:3;
% a_ell = 1:0.01:2;
for i = 1:size(a_ell,2)
    score_ell(i) = 0;
    for j = 1:size(Xagg, 2) % for each point in the point cloud,
        y_closest = fminunc(@(y) min([(max([0, Xagg(j) - -a_ell(i), -(Xagg(j) - (-3 * a_ell(i)))]))^2 + (Yagg(j) - 1)^2, ...
            (a_ell(i) * sqrt(1 - y^2) - a_ell(i) - Xagg(j))^2 + (y - Yagg(j))^2, ...
            (max([0, Xagg(j) - -a_ell(i), -(Xagg(j) - (-8 * a_ell(i)))]))^2 + (Yagg(j) - -1)^2]), Yagg(j));
        score_ell(i) = score_ell(i) + min([(max([0, Xagg(j) - -a_ell(i), -(Xagg(j) - (-3 * a_ell(i)))]))^2 + (Yagg(j) - 1)^2, ...
            (a_ell(i) * sqrt(1 - y_closest^2) - a_ell(i) - Xagg(j))^2 + (y_closest - Yagg(j))^2, ...
            (max([0, Xagg(j) - -a_ell(i), -(Xagg(j) - (-8 * a_ell(i)))]))^2 + (Yagg(j) - -1)^2]);
    end
end

fit_hyp = a_hyp(find(score_hyp == min(score_hyp)));
fit_ell = a_ell(find(score_ell == min(score_ell)));

disp(['The best hyphoid approximation is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(Xagg, 2)), ' between each point and the fitted curve.']);
disp(['The best ellipsoid approximation is with a = ', num2str(fit_ell), ', with average squared distance ', num2str(min(score_ell) / size(Xagg, 2)), ' between each point and the fitted curve.']);

hold on; plot(Xagg, Yagg, 'x');
a = fit_hyp; plot(pi / a * (-0.94:0.001:0.94) .* cot(pi * (-0.94:0.001:0.94)) - 1/a, (-0.94:0.001:0.94) / base_height(fit_hyp), 'LineWidth', 2.0); 
a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
daspect([1 1 1]);
legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(Xagg, 2))], ...
    ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(Xagg, 2))], 'Location', 'west');
exportgraphics(gcf, ['media/rojas_hyphoid_ellipsoid_fit_cell_', num2str(outline_num), '.png']);

1;

function mid = base_height(a)
    % returns the y s.t. pi*y/a * cot*pi*y) = -5. this is the base of the
    % cell
    f = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
    yl = 0; % bisection method
    yr = 0.5;
    while f(yr) > -5
        yr = (1 + yr) / 2;
    end
    mid = (yl+yr)/2;
    fmid = f(mid);
    while abs(fmid - -5) > 0.001
        if fmid < -5
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f(mid);
    end
end