load('../cell-profiles/PollenOutlines.mat');
indices = [1 5 11 17 21 26 32 37 41 46];
tipInds = [876 825 820 832 815 757 809 774 780 745]; % indices of the tip points for each outline in Outlines{1}, computed by finding a point of symmetry for curvature
% tipInds = [790 747 772 744 762 730 786 735 777 755]; % indices of the tip points for each outline in Outlines{5}, computed by finding a point of symmetry for curvature
% indices = 1:5:46; 

outline_num = 1;

hold on;
for ind = 1:size(indices, 2)
    X = Outlines_X_Lily{outline_num}(:,indices(ind));
    X = X(1:max(find(X)));
    Y = Outlines_Y_Lily{outline_num}(:,indices(ind));
    Y = Y(1:max(find(Y)));
%     plot(X, Y);
    Xtmp = X - X(tipInds(ind));
    Ytmp = Y - Y(tipInds(ind));
    
%     [ks ktheta] = compute_curvatures([X'; Y']);
%     plot(X(1:end-1)' - max(X), ks);
%     plot(ks);
%     ylim([-0.05 0.05]);
    
    angle = atan((X(tipInds(ind)-1)-X(tipInds(ind)+1))/(Y(tipInds(ind)-1)-Y(tipInds(ind)+1)));    
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    rottmp = rot * [Xtmp Ytmp]';
    Xnew = rottmp(1,:)';
    Ynew = rottmp(2,:)';
    plot(Xnew, Ynew, 'LineWidth', 2.0);

    tip_ind = find(Xnew == max(Xnew));
    rv_tmp = [Xnew(end:-1:tip_ind)' - min(Xnew(end:-1:tip_ind)); Ynew(end:-1:tip_ind)' - min(Ynew(end:-1:tip_ind))];
    rv_tmp = rv_tmp / max(rv_tmp(2,:));
    rvs{ind} = rv_tmp;
end

for ind = 1:10  % reduce to a coarser discretization of about 70 points
    rvs{ind} = rvs{ind}(:,mod(size(rvs{ind},2)-1,floor((size(rvs{ind},2)-1)/70))+1:floor((size(rvs{ind},2)-1)/70):size(rvs{ind},2));
end

for ind = 1:10
    Xagg = [Xagg rvs{ind}(1,:)];
    Yagg = [Yagg rvs{ind}(2,:)];
end

save('../cell-profiles/rojas-pollen-tube-processed.mat', 'rvs');