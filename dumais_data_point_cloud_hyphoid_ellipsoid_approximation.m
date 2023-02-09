clear all;
load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
% inds = size(Px, 2) - [0 5 10 16 20 26 30 35 40 45]; % not allunimodal
% tipX = [678.2824  645.2651  612.3227  573.9185  545.8388  505.2114  480.8379  451.4691  414.7958  383.2193];
inds = size(Px, 2) - [0 5 10 17 20 26 30 34 40 45]; % all unimodal distributions
tipX = [678.2824  645.2651  612.3227  565.7449  545.8388  505.2114  480.8379  457.1793  414.7958  383.2193];
tipY = [209.5303  208.1344  211.3621  207.7728  210.5749  201.4695  196.4848  199.1269  201.5734  198.9826];
X = [];
Y = [];
% inds = size(Px,2):-5:size(Px,2)-20;

% hold on;
for i = 1:10 % add in all the points from the five outlines
    Xnew = Px(1:max(find(Px(:,inds(i)))),inds(i))'-tipX(i);
    Ynew = Py(1:max(find(Px(:,inds(i)))),inds(i))';
%     Ynew = Py(1:max(find(Px(:,inds(i)))),inds(i))'-tipY(i);
    Indfirst = min(find(atan(diff(Ynew) ./ diff(Xnew)) >= 0)); % indices of the boundary points - exclude points where the outline curves in at the shank
    Indlast = max(find(atan(diff(Ynew) ./ diff(Xnew)) <= 0)) + 1;
    X = [X Xnew(Indfirst:Indlast)]; % the x and y coordinates
    Y = [Y Ynew(Indfirst:Indlast)]; % they trace the profile from bottom left, to tip, to top left
    
%     plot(Xnew(Indfirst:Indlast), Ynew(Indfirst:Indlast), 'x');
%     plot(Px(1:max(find(Px(:,inds(i)))),inds(i))'-tipX(i), Py(1:max(find(Px(:,inds(i)))),inds(i))', 'x');
end

Ymin = min(Y);
Ymax = max(Y);
% Ylower = Ymin * 7/8 + Ymax * 1/8;
% Yupper = Ymin * 1/8 + Ymax * 7/8;
% fitInds = find((Y >= Ylower) & (Y <= Yupper));
fitInds = find(X >= -300);
% Ylower = min(Y(fitInds));
% Yupper = max(Y(fitInds));
Ylower = min(Y);
Yupper = max(Y);

% Center around z-axis and normalize to about unit radius
X = X / (Ymax - Ymin);
Y = Y / (Ymax - Ymin);
Y = Y - mean([Ymax, Ymin]) / (Ymax - Ymin);
X = 2 * X;
Y = 2 * Y;

% % ===== TEST =====
% % use only the bottom half of the outline
% X = X(find(Y >= 0));
% Y = Y(find(Y >= 0));
% % ===== END TEST =====

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
a_ell = 1.4:0.01:2.4;
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

disp(['The best hyphoid approximation is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(X,2)), ' between each point and the fitted curve.']);
disp(['The best ellipsoid approximation is with a = ', num2str(fit_ell), ', with average squared distance ', num2str(min(score_ell) / size(X,2)), ' between each point and the fitted curve.']);

hold on; plot(X, Y, 'x');
a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96) / base_height(fit_hyp), 'LineWidth', 2.0); 
% a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
set(gca, 'ColorOrderIndex', 3);
plot([pi / fit_hyp * 0.96 * cot(pi * 0.96) - 1/fit_hyp, -fit_ell], [1, 1], 'LineWidth', 2.0); % cylinder region for the ellipsoid tip
set(gca, 'ColorOrderIndex', 3);
plot([pi / fit_hyp * 0.96 * cot(pi * 0.96) - 1/fit_hyp, -fit_ell], [-1, -1], 'LineWidth', 2.0);
daspect([1 1 1]);
legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(X, 2))], ...
    ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(X, 2))], 'Location', 'west');
ax = gca;
set(ax, 'fontsize', 12);
xlim tight; ylim([-1.1 1.1]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4]);
pause(0.5);
exportgraphics(gcf, ['media/dumais_hyphoid_ellipsoid_fit.png']);
close all;

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

%% Old polynomial fitting

% % degree = 4
% % Xp = polyfit(Y(fitInds), X(fitInds), 4);
% % Xpoly = @(t) Xp(1) * t.^4 + Xp(2) * t.^3 + Xp(3) * t.^2+ Xp(4) * t + Xp(5);
% % degree = 6
% % Xp = polyfit(Y(fitInds), X(fitInds), 6);
% % Xpoly = @(t) Xp(1) * t.^6 + Xp(2) * t.^5 + Xp(3) * t.^4+ Xp(4) * t.^3 + Xp(5) * t.^2 + Xp(6) * t + Xp(7);
% % degree = 8
% % Xp = polyfit(Y(fitInds), X(fitInds), 8);
% % Xpoly = @(t) Xp(1) * t.^8 + Xp(2) * t.^7 + Xp(3) * t.^6 + Xp(4) * t.^5 + Xp(5) * t.^4 + Xp(6) * t.^3 + Xp(7) * t.^2 + Xp(8) * t + Xp(9);
% % degree = 10
% Xp = polyfit(Y(fitInds), X(fitInds), 10);
% Xpoly = @(t) Xp(1) * t.^10 + Xp(2) * t.^9 + Xp(3) * t.^8 + Xp(4) * t.^7 + Xp(5) * t.^6 + Xp(6) * t.^5 + Xp(7) * t.^4 + Xp(8) * t.^3 + Xp(9) * t.^2 + Xp(10) * t + Xp(11);
% % degree = 12;
% % Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 12);
% % Xpoly = @(t) Xp(1) * t.^12 + Xp(2) * t.^11 + Xp(3) * t.^10 + Xp(4) * t.^9 + Xp(5) * t.^8 + Xp(6) * t.^7 + Xp(7) * t.^6 + Xp(8) * t.^5 + Xp(9) * t.^4 + Xp(10) * t.^3 + Xp(11) * t.^2 + Xp(12) * t + Xp(13);
% % degree = 14;
% % Xp = polyfit(Y(floor(size(X,2))/6:5*floor(size(X,2))/6), X(floor(size(X,2))/6:5*floor(size(X,2))/6), 14);
% % Xpoly = @(t) Xp(1) * t.^14 + Xp(2) * t.^13 + Xp(3) * t.^12 + Xp(4) * t.^11 + Xp(5) * t.^10 + Xp(6) * t.^9 + Xp(7) * t.^8 + Xp(8) * t.^7 + Xp(9) * t.^6 + Xp(10) * t.^5 + Xp(11) * t.^4 + Xp(12) * t.^3 + Xp(13) * t.^2 + Xp(14) * t + Xp(15);
% 
% yrng = Ylower:1:Yupper;
% for i = 1:size(yrng,2)
%     % degree = 8
% %     if yrng(i) < mean([Ylower Yupper])
% %         s(i) = -midpt(@(t) sqrt(1 + (8 * Xp(1) * t.^7 + 7 * Xp(2) * t.^6 + 6 * Xp(3) * t.^5 + 5 * Xp(4) * t.^4 + 4 * Xp(5) * t.^3 + 3 * Xp(6) * t.^2 + 2 * Xp(7) * t + Xp(8))^2), yrng(i), mean([Ylower Yupper]), 256);
% %     else
% %         s(i) = midpt(@(t) sqrt(1 + (8 * Xp(1) * t.^7 + 7 * Xp(2) * t.^6 + 6 * Xp(3) * t.^5 + 5 * Xp(4) * t.^4 + 4 * Xp(5) * t.^3 + 3 * Xp(6) * t.^2 + 2 * Xp(7) * t + Xp(8))^2), mean([Ylower Yupper]), yrng(i), 256);
% %     end
% %     y = yrng(i);
% %     ks(i) = - sum([8:-1:2] .* [7:-1:1] .* Xp(1:end-2) .* y .^ ([6:-1:0])) / (1 + (sum([8:-1:1] .* Xp(1:end-1) .* y .^ (7:-1:0)))^2)^(3/2);
%     % degree = 10
%     if yrng(i) < mean([Ylower Yupper])
%         s(i) = -midpt(@(t) sqrt(1 + (10 * Xp(1) * t.^9 + 9 * Xp(2) * t.^8 + 8 * Xp(3) * t.^7 + 7 * Xp(4) * t.^6 + 6 * Xp(5) * t.^5 + 5 * Xp(6) * t.^4 + 4 * Xp(7) * t.^3 + 3 * Xp(8) * t.^2 + 2 * Xp(9) * t + Xp(10))^2), yrng(i), mean([Ylower Yupper]), 256);
%     else
%         s(i) = midpt(@(t) sqrt(1 + (10 * Xp(1) * t.^9 + 9 * Xp(2) * t.^8 + 8 * Xp(3) * t.^7 + 7 * Xp(4) * t.^6 + 6 * Xp(5) * t.^5 + 5 * Xp(6) * t.^4 + 4 * Xp(7) * t.^3 + 3 * Xp(8) * t.^2 + 2 * Xp(9) * t + Xp(10))^2), mean([Ylower Yupper]), yrng(i), 256);
%     end
%     y = yrng(i);
%     ks(i) = - sum([10:-1:2] .* [9:-1:1] .* Xp(1:end-2) .* y .^ ([8:-1:0])) / (1 + (sum([10:-1:1] .* Xp(1:end-1) .* y .^ (9:-1:0)))^2)^(3/2);
% end
% 
% hold on; 
% plot(Xpoly(Ylower:0.01:Yupper), Ylower:0.01:Yupper);
% plot(X, Y, 'x');
% daspect([1 1 1]);
% close all;