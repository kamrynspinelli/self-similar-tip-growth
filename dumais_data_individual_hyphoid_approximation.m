clear all;
load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
inds = size(Px, 2) - [0 5 10 17 20 26 30 34 40 45]; % all unimodal distributions
tipX = [678.2824  645.2651  612.3227  565.7449  545.8388  505.2114  480.8379  457.1793  414.7958  383.2193];

for prof = 1:10
    % get one outline
    X = Px(1:max(find(Px(:,inds(prof)))),inds(prof))'-tipX(prof); % the x and y coordinates
    Y = Py(1:max(find(Px(:,inds(prof)))),inds(prof))'; % they trace the profile from bottom left, to tip, to top left
    
    Indfirst = min(find(atan(diff(Y) ./ diff(X)) >= 0)); % indices of the boundary points - exclude points where the outline curves in at the shank
    Indlast = max(find(atan(diff(Y) ./ diff(X)) <= 0)) + 1;
    X = X(Indfirst:Indlast);
    Y = Y(Indfirst:Indlast);
    
    Ymin = min(Y); % tip-align and normalize all the points
    Ymax = max(Y);
%     fitInds = find(X >= -300);
%     Ylower = min(Y(fitInds));
%     Yupper = max(Y(fitInds));
    Ylower = min(Y);
    Yupper = max(Y);
    X = X / (Ymax - Ymin);
    Y = Y / (Ymax - Ymin);
    Y = Y - mean([Ymax, Ymin]) / (Ymax - Ymin);
    X = 2 * X;
    Y = 2 * Y;
    
%     % ===== TEST =====
%     % use only the bottom half of the outline
%     X = X(find(Y >= 0));
%     Y = Y(find(Y >= 0));
%     % ===== END TEST =====
    
    % do the hyphoid fitting
    a_hyp = 2:0.1:12;
    for i = 1:size(a_hyp,2)
        score_hyp(i) = 0;
        for j = 1:size(X, 2) % for each point in the point cloud,
%             y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y - Y(j))^2, 0.0001);
%             score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y_closest - Y(j))^2;
            y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (y - Y(j))^2, 0.0001);
            score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (y_closest - Y(j))^2;
        end
    end
    fit_hyp(prof) = a_hyp(find(score_hyp == min(score_hyp)));
end

save('dumais-data-individual-hyphoid-fits.mat', 'fit_hyp');

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