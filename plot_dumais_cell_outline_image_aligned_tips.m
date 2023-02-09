load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
hold on;
c = lines(8); % some colors to use
for profile_num = 1:5
    profile = size(Px,2) - 0 - 5 * (profile_num - 1);
    X = Px(1:max(find(Px(:,profile))),profile)'; % load the x and y coordinates
    Y = Py(1:max(find(Px(:,profile))),profile)'; % they trace the profile from bottom left, to tip, to top left
    Y = -Py(1:max(find(Px(:,profile))),profile)'; % use the lower profile instead of the top
    Y = Y - min(Y) + 1;
    [AX BX CX DX] = find_splines(1:size(X,2), X, 1, -1); % interpolate the x and y coordinates with cubic splines
    [AY BY CY DY] = find_splines(1:size(X,2), Y, 0, 0);
    
    i_maxX = find(X == max(X)); % find the index of the rightmost point
    if X(i_maxX-1) > X(i_maxX+1) % find the index of the next-rightmost point
        i_2ndmaxX = i_maxX-1;
    else
        i_2ndmaxX = i_maxX+1;
    end
    % find the 'index' of the interpolated tip using quadratic formula
    if abs(i_maxX - (-2 * BX(min(i_maxX, i_2ndmaxX)) - sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)))) < 1
        i_tip = (-2 * BX(min(i_maxX, i_2ndmaxX)) - sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)));
    else
        i_tip = (-2 * BX(min(i_maxX, i_2ndmaxX)) + sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)));
    end
    tipX(profile_num) = AX(floor(i_tip)) * i_tip^3 + BX(floor(i_tip)) * i_tip^2 + CX(floor(i_tip)) * i_tip + DX(floor(i_tip));
    tipY(profile_num) = AY(floor(i_tip)) * i_tip^3 + BY(floor(i_tip)) * i_tip^2 + CY(floor(i_tip)) * i_tip + DY(floor(i_tip));
    
    for i = 1:size(X,2)-1
        splX = @(x) AX(i) * x .^3 + BX(i) * x .^2 + CX(i) * x + DX(i);
        splY = @(x) AY(i) * x .^3 + BY(i) * x .^2 + CY(i) * x + DY(i);
        p = plot(splX(i:0.001:i+1) - tipX(profile_num), splY(i:0.001:i+1) - tipY(profile_num), 'LineWidth', 2.0);
        p.Color = c(profile_num,:);
    end
end
daspect([1 1 1]);
set(gca,'visible','off');
exportgraphics(gcf, 'media/dumais-root-hair-profiles-aligned-tips.png');