clear all;
load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
% inds = size(Px, 2) - [0 5 10 16 20 26 30 35 40 45]; % not allunimodal
% tipX = [678.2824  645.2651  612.3227  573.9185  545.8388  505.2114  480.8379  451.4691  414.7958  383.2193];
inds = size(Px, 2) - [0 5 10 17 20 26 30 34 40 45]; % all unimodal distributions
tipX = [678.2824  645.2651  612.3227  565.7449  545.8388  505.2114  480.8379  457.1793  414.7958  383.2193];
X = [];
Y = [];
% inds = size(Px,2):-5:size(Px,2)-20;

% hold on;
for i = 1:10 % add in all the points from the five outlines
    Xnew = Px(1:max(find(Px(:,inds(i)))),inds(i))'-tipX(i);
    Ynew = Py(1:max(find(Px(:,inds(i)))),inds(i))';
    
    Indfirst = min(find(atan(diff(Ynew) ./ diff(Xnew)) >= 0)); % indices of the boundary points - exclude points where the outline curves in at the shank
    Indlast = max(find(atan(diff(Ynew) ./ diff(Xnew)) <= 0)) + 1;
    Xnew = Xnew(Indfirst:Indlast);
    Ynew = Ynew(Indfirst:Indlast);

    Ymin = min(Ynew);
    Ymax = max(Ynew);

    % Center around z-axis and normalize to about unit radius
    Xnew = Xnew / (Ymax - Ymin);
    Ynew = Ynew / (Ymax - Ymin);
    Ynew = Ynew - mean([Ymax, Ymin]) / (Ymax - Ymin);
    Xnew = 2 * Xnew;
    Ynew = 2 * Ynew;
    X = [X, Xnew]; % the x and y coordinates
    Y = [Y, Ynew]; % they trace the profile from bottom left, to tip, to top left
    
%     plot(Px(1:max(find(Px(:,inds(i)))),inds(i))'-tipX(i), Py(1:max(find(Px(:,inds(i)))),inds(i))');
end

Ymin = min(Y);
Ymax = max(Y);
% Ylower = Ymin * 7/8 + Ymax * 1/8;
% Yupper = Ymin * 1/8 + Ymax * 7/8;
% fitInds = find((Y >= Ylower) & (Y <= Yupper));
fitInds = find(X >= -300);
Ylower = min(Y(fitInds));
Yupper = max(Y(fitInds));

% find the hyphoid and ellipsoid shapes minimizing the sum of squared
% distances to the points
a_hyp = 2:0.1:8;
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

% a_ell = 0.5:0.1:3;
% a_ell = 1.4:0.01:2.4;
% for i = 1:size(a_ell,2)
%     score_ell(i) = 0;
%     for j = 1:size(X, 2) % for each point in the point cloud,
%         y_closest = fminunc(@(y) min([(max([0, X(j) - -a_ell(i), -(X(j) - (-3 * a_ell(i)))]))^2 + (Y(j) - 1)^2, ...
%             (a_ell(i) * sqrt(1 - y^2) - a_ell(i) - X(j))^2 + (y - Y(j))^2, ...
%             (max([0, X(j) - -a_ell(i), -(X(j) - (-8 * a_ell(i)))]))^2 + (Y(j) - -1)^2]), Y(j));
%         score_ell(i) = score_ell(i) + min([(max([0, X(j) - -a_ell(i), -(X(j) - (-3 * a_ell(i)))]))^2 + (Y(j) - 1)^2, ...
%             (a_ell(i) * sqrt(1 - y_closest^2) - a_ell(i) - X(j))^2 + (y_closest - Y(j))^2, ...
%             (max([0, X(j) - -a_ell(i), -(X(j) - (-8 * a_ell(i)))]))^2 + (Y(j) - -1)^2]);
%     end
% end

fit_hyp = a_hyp(find(score_hyp == min(score_hyp)));
% fit_ell = a_ell(find(score_ell == min(score_ell)));

disp(['The best hyphoid approximation is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(X,2)), ' between each point and the fitted curve.']);
% disp(['The best ellipsoid approximation is with a = ', num2str(fit_ell), ', with average squared distance ', num2str(min(score_ell) / size(X,2)), ' between each point and the fitted curve.']);

hold on; plot(X, Y, 'x');
% a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96) / base_height(fit_hyp), 'LineWidth', 2.0); 
a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96), 'LineWidth', 2.0); 
% a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
% daspect([1 1 1]);
% legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(X, 2))], ...
%     ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(X, 2))], 'Location', 'west');
% exportgraphics(gcf, ['media/dumais_hyphoid_ellipsoid_fit.png']);
% close all;

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