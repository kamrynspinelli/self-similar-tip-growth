clear all;
load('../cell-profiles/PollenOutlines.mat');
indices = [1 5 11 17 21 26 32 37 41 46];
% indices = 1:5:46;
% tipInds = [876 825 820 832 815 757 809 774 780 745]; % indices of the tip points for each outline, computed by finding a point of symmetry for curvature

outline_num = 1;

kssum = zeros(1,60);
kthetasum = zeros(1,60);

Xagg = [];
Yagg = [];

hold on;
for ind = 1:size(indices, 2)
    X = Outlines_X_Lily{outline_num}(:,indices(ind));
    X = X(1:max(find(X)));
    Y = Outlines_Y_Lily{outline_num}(:,indices(ind));
    Y = Y(1:max(find(Y)));
    X = X - min(X);
    Y = Y - min(Y);
    
%     ext = 2000;
%     augX(ext+1:ext+size(X,1)) = X;
%     augY(ext+1:ext+size(Y,1)) = Y;
%     for i = 1:ext
%         augX(ext+1-i) = augX(ext+1-i+1) - (X(2) - X(1));
%         augY(ext+1-i) = augY(ext+1-i+1) - (Y(2) - Y(1));
%         augX(ext+size(X,1)+i) = augX(ext+size(X,1)+i-1) + (X(end) - X(end-1));
%         augY(ext+size(X,1)+i) = augY(ext+size(X,1)+i-1) + (Y(end) - Y(end-1));
%     end

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
    a = -1/m;
    c = yint/m - xint;
%     close all; hold on; plot(X, Y); plot(X(tip_ind), Y(tip_ind), 'x'); daspect([1 1 1])
%     plot([0 300], m * ([0 300] - X(tip_ind)) + Y(tip_ind));
    
%     for i = 1:size(X,1)
%         a = - (mean(augX) - augX(i)) / (mean(augY) - augY(i));
%         c = -mean(augX + a * augY);
%         score(i) = (a^2 + 1)^-1 * sum((augX + a * augY + c).^2);
%     end
    
%     coords = [X(:), Y(:)];
%     [U, S, V] = svd(coords - repmat(mean(coords), size(coords, 1), 1), 0);
%     a = -V(1, end) / V(2, end);
%     b = mean(coords * V(:, end)) / V(2, end);
    
%     tip_ind = find(diff(score(1:end-1)) .* diff(score(2:end)) < 0);
%     tip_ind = tip_ind(2) + 1; % get the tip point
%     tip_ind = find(score == min(score));
%     a = - (mean(augX) - augX(tip_ind)) / (mean(augY) - augY(tip_ind));
%     c = -mean(augX + a * augY);
    angle = - atan(-1/a);
%     angle = atan((X(tip_ind-1)-X(tip_ind+1))/(Y(tip_ind-1)-Y(tip_ind+1)));
    Xtmp = X - X(tip_ind);
    Ytmp = Y - Y(tip_ind);
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    rottmp = rot * [Xtmp Ytmp]';
    Xnew = rottmp(1,:)';
    Ynew = rottmp(2,:)';
    plot(Xnew, Ynew);
    
    tip_ind = find(Xnew == max(Xnew));
    Xnew = Xnew - Xnew(tip_ind);
    Ynew = Ynew - Ynew(tip_ind);
%     hold on; plot(Xnew, Ynew, 'LineWidth', 2.0);
%     plot(-300:0, zeros(1,301)); daspect([1 1 1]);
    
    rv_tmp = [Xnew(end:-1:tip_ind)' - min(Xnew(end:-1:tip_ind)); Ynew(end:-1:tip_ind)' - min(Ynew(end:-1:tip_ind))];
    rest_tmp = [Xnew(1:tip_ind)' - min(Xnew(end:-1:tip_ind)); -(Ynew(1:tip_ind)' - min(Ynew(end:-1:tip_ind)))];
    rv_tmp = rv_tmp / max(rv_tmp(2,:));
    rest_tmp = rest_tmp / max(rest_tmp(2,:));
    rvs{ind} = rv_tmp;

%     [ks1 ktheta1] = compute_curvatures_spline(rv_tmp);
%     [ks2 ktheta2] = compute_curvatures_spline(rest_tmp);
%     
%     plot(-size(ks2,2):size(ks1,2)-1, [ks2 fliplr(ks1)]); % all the points
%     plot(-size(ktheta2,2):size(ktheta1,2)-1, [ktheta2 fliplr(ktheta1)]); % all the points

    [ks1 ktheta1] = compute_curvatures(rv_tmp(:,mod(size(rv_tmp,2)-1,floor((size(rv_tmp,2)-1)/30))+1:floor((size(rv_tmp,2)-1)/30):size(rv_tmp,2)));
    [ks2 ktheta2] = compute_curvatures(rest_tmp(:,mod(size(rest_tmp,2)-1,floor((size(rest_tmp,2)-1)/30))+1:floor((size(rest_tmp,2)-1)/30):size(rest_tmp,2)));
    
%     close all; hold on;
%     plot(-size(ks2,2):size(ks1,2)-1, [ks2 fliplr(ks1)]); % trimmed subset
%     plot(-size(ktheta2,2):size(ktheta1,2)-1, [ktheta2 fliplr(ktheta1)]); % trimmed subset
    
%     close all; hold on;
%     plot(rv_tmp(1,:) - max(rv_tmp(1,:)), rv_tmp(2,:), rest_tmp(1,:) - max(rest_tmp(1,:)), rest_tmp(2,:)); daspect([1 1 1]);

    kssum = kssum + [ks2(1:30) fliplr(ks1(end-29:end))];
    kthetasum = kthetasum + [ktheta2(1:30) fliplr(ktheta1(end-29:end))];
end

for ind = 1:10  % reduce to a coarser discretization of about 70 points
    rvs{ind} = rvs{ind}(:,mod(size(rvs{ind},2)-1,floor((size(rvs{ind},2)-1)/70))+1:floor((size(rvs{ind},2)-1)/70):size(rvs{ind},2));
end



for ind = 1:10
    Xagg = [Xagg rvs{ind}(1,:)];
    Yagg = [Yagg rvs{ind}(2,:)];
end

save('../cell-profiles/rojas-pollen-tube-processed.mat', 'rvs');