clear all;

N = 128; % number of points to put on the discretized outline

%% Caulonema
% convert the text file to coordinates
filename = '../cell-profiles/moss/Caulonema control/Snakes-Caulonema-Control-trimmed.txt'; 
text = fileread(filename);
splitted = split(text, '#');
for i=1:size(splitted,1)
    scanned = textscan(splitted{i}, '%f');
    coords{i}(1,:) = scanned{1}(4:5:end)';
    coords{i}(2,:) = scanned{1}(5:5:end)';
%     plot(coords{i}(1,:)-max(coords{i}(1,:)), coords{i}(2,:)-coords{i}(2,find(coords{i}(1,:) == max(coords{i}(1,:))))); daspect([1 1 1]);
end
% select some of the outlines that look more or less symmetric
% inds = [1 2 4 5 7 14 15 22];
inds = [1 2 5 7 14 15];
for i = 1:size(inds, 2)
%     profiles{i} = coords{inds(i)}(:,find(coords{inds(i)}(1,:) >= 250));
    profiles{i} = coords{inds(i)}(:, ...
        find(abs(coords{inds(i)}(1,:) - max(coords{inds(i)}(1,:))) <= 120));
end

% rotate the outlines
for i = 1:size(inds,2)
    clear rot_smple_points;
    rot_smple_points = zeros(2,1);
    spot = 1; % keep track of how many sample points we're using for the least squares computation
    for j = 1:size(profiles{i},2) % build the list of points on the tube and old end
        rot_smple_points(1,spot) = profiles{i}(1,j);
        rot_smple_points(2,spot) = profiles{i}(2,j);
        spot = spot + 1;
    end
    coords = [rot_smple_points(1,:)', rot_smple_points(2,:)']; % SVD magic to find the line minimizing the sum of squared distances: https://stackoverflow.com/a/12251893
    [U, S, V] = svd(coords - repmat(mean(coords), size(coords, 1), 1), 0);
    a = -V(1, end) / V(2, end);
    b = mean(coords * V(:, end)) / V(2, end);
    angle = - atan(-1/a) + pi/2;
    if angle > pi/2
        angle = angle - pi;
    end
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    rottmp = rot * [profiles{i}(1,:)' profiles{i}(2,:)'-b]'; % [Cx(:,i) Cy(:,i)-b]'
    profiles_rot{i}(1,:) = rottmp(1,:);
    profiles_rot{i}(2,:) = rottmp(2,:);
%     plot(profiles_rot{i}(1,:)-max(profiles_rot{i}(1,:)), profiles_rot{i}(2,:)-profiles_rot{i}(2,find(profiles_rot{i}(1,:) == max(profiles_rot{i}(1,:))))); daspect([1 1 1])
end

% hold on; for i = 1:size(inds,2); plot(profiles_rot{i}(1,:)-max(profiles_rot{i}(1,:)), profiles_rot{i}(2,:)-profiles_rot{i}(2,find(profiles_rot{i}(1,:) == max(profiles_rot{i}(1,:))))); end; daspect([1 1 1]);
% hold on; for i = 1:size(inds,2); p = profiles_rot{i}; [ks ktheta] = compute_curvatures(p); plot(cumsum(sqrt(diff(p(2,1:end-1))).^2 + (diff(p(1,1:end-1))).^2), ks(1:end-1)); end;

for prof = 1:size(inds,2)
% for prof = 1:size(Cx,1)
    X = profiles_rot{prof}(1,:); % load the x and y coordinates for this profile. 1:32 selects the old end of the yeast
    Y = profiles_rot{prof}(2,:); % they trace the profile from bottom left, to tip, to top left
%     Y = -Y; % use the lower profile instead of the top
%     X = fliplr(X);
%     Y = Y - min(Y) + 1;
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
    % compute the total arclength from the tip to the top left endpoint
    total_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
        + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
    for i = ceil(i_tip):size(X,2)-1 % add up the arclengths for all the patches from the tip to the rear
        total_arclength = total_arclength ...
            + midpt(@(t) sqrt((3 * AX(i) * t.^2 + 2 * BX(i) * t + CX(i))^2 ...
                + (3 * AY(i) * t.^2 + 2 * BY(i) * t + CY(i))^2), i, i+1, 64);
    end
    target_arclengths = fliplr(0:total_arclength/N:total_arclength); % the desired arclength from the tip to the i-th marker point

    rv(1,1) = X(end); rv(2,1) = Y(end); % preload the tip and rear points
    rv(1,N+1) = AX(floor(i_tip)) * i_tip^3 + BX(floor(i_tip)) * i_tip^2 + CX(floor(i_tip)) * i_tip + DX(floor(i_tip));
    rv(2,N+1) = AY(floor(i_tip)) * i_tip^3 + BY(floor(i_tip)) * i_tip^2 + CY(floor(i_tip)) * i_tip + DY(floor(i_tip));
    for i = 2:N % find the positions of the rest of the points
        targ = target_arclengths(i);
        arclength_sum = 0;
        interp = floor(i_tip); % the number of the tipmost spline
        tip_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
        + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
        % find which patch has the desired arclength coordinate
        if targ > tip_arclength
            % if the target arclength is less than the arclength from the tip
            % to the adjacent marker point, need to advance through the other
            % patches
            arclength_sum = arclength_sum + tip_arclength;
            interp = floor(i_tip) + 1;
            while arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64) < targ
                arclength_sum = arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64);
                interp = interp + 1;
            end
            remaining_arclength = targ - arclength_sum;
            % use bisection method to find the point on this patch which gives the
            % desired arclength
            i_l = interp;
            i_r = interp + 1;
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
            while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
                if guess_arclength - remaining_arclength < 0
                    i_l = mid;
                else
                    i_r = mid;
                end
                mid = (i_l + i_r) / 2;
                guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
            end
            rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
            rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
        else % if this arclength coordinate lies on the tip interpolant, just do bisection method
            remaining_arclength = targ;
            % use bisection method to find the point on this patch which gives the
            % desired arclength
            i_l = i_tip;
            i_r = ceil(i_tip);
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
            while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
                if guess_arclength - remaining_arclength < 0
                    i_l = mid;
                else
                    i_r = mid;
                end
                mid = (i_l + i_r) / 2;
                guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
            end
            rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
            rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
        end
    end
%     rv(2,:) = -rv(2,:); % use the lower profile instead of the top
    rv(2,:) = rv(2,:) - min(rv(2,:)); % center at (0,0)
    rv(1,:) = rv(1,:) - min(rv(1,:));
    rv = rv / max(rv(2,:)); % normalize to have height 1
    rvs{prof} = rv;
%     hold on; plot(rv(1,:)-rv(1,end), rv(2,:)); daspect([1 1 1]);
    
%     [ks ktheta] = compute_curvatures_spline(rv);
%     ksuppers{prof} = ks;
%     kthetauppers{prof} = ktheta;
%     
%     kssumupper = kssumupper + ks;
%     kthetasumupper = kthetasumupper + ktheta;
%     
%     plot(0:N, fliplr(ks), 0:N, fliplr(ktheta));
end

save('../cell-profiles/moss/caulonema-data.mat', 'rvs');

% X = [];
% Y = [];
% for i = 1:size(inds, 2)
%     X = [X rvs{i}(1,1:4:end)-max(rvs{i}(1,1:4:end))];
%     Y = [Y rvs{i}(2,1:4:end)];
% end
% 
% a_hyp = 2:0.1:8;
% % a_hyp = 5:0.01:6.5;
% for i = 1:size(a_hyp,2)
%     score_hyp(i) = 0;
%     for j = 1:size(X, 2) % for each point in the point cloud,
% %         y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y - Y(j))^2, 0.0001);
% %         score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (1 / base_height(a_hyp(i)) * y_closest - Y(j))^2;
%         y_closest = fminunc(@(y) (pi * y / a_hyp(i) * cot(pi * y) - 1/a_hyp(i) - X(j))^2 + (y - Y(j))^2, 0.0001);
%         score_hyp(i) = score_hyp(i) + (pi * y_closest / a_hyp(i) * cot(pi * y_closest) - 1/a_hyp(i) - X(j))^2 + (y_closest - Y(j))^2;
%     end
% end
% 
% fit_hyp = a_hyp(find(score_hyp == min(score_hyp)));
% 
% disp(['The best hyphoid approximation for caulonema is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(X,2)), ' between each point and the fitted curve.']);
% 
% hold on; plot(X, Y, 'x');
% a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96) / base_height(fit_hyp), 'LineWidth', 2.0); 
% % a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
% daspect([1 1 1]);
% legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(X, 2))], 'Location', 'west');
% %     ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(X, 2))], 'Location', 'west');
% ax = gca;
% set(ax, 'fontsize', 12);
% exportgraphics(gcf, ['media/moss_caulonema_hyphoid_fit.png']);
% close all;

%% Chloronema

clearvars -except N;

% convert the text file to coordinates
filename = '../cell-profiles/moss/Chloronema control/Snakes-Chloronema_control-trimmed.txt'; 
text = fileread(filename);
splitted = split(text, '#');
for i=1:size(splitted,1)
    scanned = textscan(splitted{i}, '%f');
    coords{i}(1,:) = scanned{1}(4:5:end)';
    coords{i}(2,:) = scanned{1}(5:5:end)';
%     plot(coords{i}(1,:)-max(coords{i}(1,:)), coords{i}(2,:)-coords{i}(2,find(coords{i}(1,:) == max(coords{i}(1,:))))); daspect([1 1 1]);
end
% select some of the outlines that look more or less symmetric
inds = [1 2 4 8 9 13 25];
for i = 1:size(inds, 2)
%     profiles{i} = coords{inds(i)}(:,find(coords{inds(i)}(1,:) >= 250));
    profiles{i} = coords{inds(i)}(:, ...
        find(abs(coords{inds(i)}(1,:) - max(coords{inds(i)}(1,:))) <= 120));
end

% rotate the outlines
for i = 1:size(inds,2)
    clear rot_smple_points;
    rot_smple_points = zeros(2,1);
    spot = 1; % keep track of how many sample points we're using for the least squares computation
    for j = 1:size(profiles{i},2) % build the list of points on the tube and old end
        rot_smple_points(1,spot) = profiles{i}(1,j);
        rot_smple_points(2,spot) = profiles{i}(2,j);
        spot = spot + 1;
    end
    coords = [rot_smple_points(1,:)', rot_smple_points(2,:)']; % SVD magic to find the line minimizing the sum of squared distances: https://stackoverflow.com/a/12251893
    [U, S, V] = svd(coords - repmat(mean(coords), size(coords, 1), 1), 0);
    a = -V(1, end) / V(2, end);
    b = mean(coords * V(:, end)) / V(2, end);
    angle = - atan(-1/a) + pi/2;
    if angle > pi/2
        angle = angle - pi;
    end
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    rottmp = rot * [profiles{i}(1,:)' profiles{i}(2,:)'-b]'; % [Cx(:,i) Cy(:,i)-b]'
    profiles_rot{i}(1,:) = rottmp(1,:);
    profiles_rot{i}(2,:) = rottmp(2,:);
end

hold on; for i = 1:size(inds,2); plot(profiles_rot{i}(1,:)-max(profiles_rot{i}(1,:)), profiles_rot{i}(2,:)-profiles_rot{i}(2,find(profiles_rot{i}(1,:) == max(profiles_rot{i}(1,:))))); end; daspect([1 1 1]);
% hold on; for i = 1:size(inds,2); p = profiles_rot{i}; [ks ktheta] = compute_curvatures(p); plot(cumsum(sqrt(diff(p(2,1:end-1))).^2 + (diff(p(1,1:end-1))).^2), ks(1:end-1)); end;

for prof = 1:size(inds,2)
% for prof = 1:size(Cx,1)
    X = profiles_rot{prof}(1,:); % load the x and y coordinates for this profile. 1:32 selects the old end of the yeast
    Y = profiles_rot{prof}(2,:); % they trace the profile from bottom left, to tip, to top left
%     Y = -Y; % use the lower profile instead of the top
%     X = fliplr(X);
%     Y = Y - min(Y) + 1;
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
    % compute the total arclength from the tip to the top left endpoint
    total_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
        + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
    for i = ceil(i_tip):size(X,2)-1 % add up the arclengths for all the patches from the tip to the rear
        total_arclength = total_arclength ...
            + midpt(@(t) sqrt((3 * AX(i) * t.^2 + 2 * BX(i) * t + CX(i))^2 ...
                + (3 * AY(i) * t.^2 + 2 * BY(i) * t + CY(i))^2), i, i+1, 64);
    end
    target_arclengths = fliplr(0:total_arclength/N:total_arclength); % the desired arclength from the tip to the i-th marker point

    rv(1,1) = X(end); rv(2,1) = Y(end); % preload the tip and rear points
    rv(1,N+1) = AX(floor(i_tip)) * i_tip^3 + BX(floor(i_tip)) * i_tip^2 + CX(floor(i_tip)) * i_tip + DX(floor(i_tip));
    rv(2,N+1) = AY(floor(i_tip)) * i_tip^3 + BY(floor(i_tip)) * i_tip^2 + CY(floor(i_tip)) * i_tip + DY(floor(i_tip));
    for i = 2:N % find the positions of the rest of the points
        targ = target_arclengths(i);
        arclength_sum = 0;
        interp = floor(i_tip); % the number of the tipmost spline
        tip_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
        + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
        % find which patch has the desired arclength coordinate
        if targ > tip_arclength
            % if the target arclength is less than the arclength from the tip
            % to the adjacent marker point, need to advance through the other
            % patches
            arclength_sum = arclength_sum + tip_arclength;
            interp = floor(i_tip) + 1;
            while arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64) < targ
                arclength_sum = arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64);
                interp = interp + 1;
            end
            remaining_arclength = targ - arclength_sum;
            % use bisection method to find the point on this patch which gives the
            % desired arclength
            i_l = interp;
            i_r = interp + 1;
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
            while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
                if guess_arclength - remaining_arclength < 0
                    i_l = mid;
                else
                    i_r = mid;
                end
                mid = (i_l + i_r) / 2;
                guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
            end
            rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
            rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
        else % if this arclength coordinate lies on the tip interpolant, just do bisection method
            remaining_arclength = targ;
            % use bisection method to find the point on this patch which gives the
            % desired arclength
            i_l = i_tip;
            i_r = ceil(i_tip);
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
            while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
                if guess_arclength - remaining_arclength < 0
                    i_l = mid;
                else
                    i_r = mid;
                end
                mid = (i_l + i_r) / 2;
                guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                    + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
            end
            rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
            rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
        end
    end
%     rv(2,:) = -rv(2,:); % use the lower profile instead of the top
    rv(2,:) = rv(2,:) - min(rv(2,:)); % center at (0,0)
    rv(1,:) = rv(1,:) - min(rv(1,:));
    rv = rv / max(rv(2,:)); % normalize to have height 1
    rvs{prof} = rv;
%     hold on; plot(rv(1,:)-rv(1,end), rv(2,:)); daspect([1 1 1]);
    
%     [ks ktheta] = compute_curvatures_spline(rv);
%     ksuppers{prof} = ks;
%     kthetauppers{prof} = ktheta;
%     
%     kssumupper = kssumupper + ks;
%     kthetasumupper = kthetasumupper + ktheta;
%     
%     plot(0:N, fliplr(ks), 0:N, fliplr(ktheta));
end

save('../cell-profiles/moss/chloronema-data.mat', 'rvs');

X = [];
Y = [];
for i = 1:size(inds, 2)
    X = [X rvs{i}(1,1:4:end)-max(rvs{i}(1,1:4:end))];
    Y = [Y rvs{i}(2,1:4:end)];
end

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

fit_hyp = a_hyp(find(score_hyp == min(score_hyp)));

disp(['The best hyphoid approximation for chloronema is with a = ', num2str(fit_hyp), ', with average squared distance ', num2str(min(score_hyp) / size(X,2)), ' between each point and the fitted curve.']);

hold on; plot(X, Y, 'x');
a = fit_hyp; plot(pi / a * (-0.96:0.001:0.96) .* cot(pi * (-0.96:0.001:0.96)) - 1/a, (-0.96:0.001:0.96) / base_height(fit_hyp), 'LineWidth', 2.0); 
% a = fit_ell; plot(a * sqrt(1 - (-1:0.001:1).^2) - a, (-1:0.001:1), 'LineWidth', 2.0); daspect([1 1 1]);
daspect([1 1 1]);
legend('Point cloud', ['Hyphoid fit: a = ', num2str(fit_hyp), ', avg. sq. dist. = ', num2str(min(score_hyp) / size(X, 2))], 'Location', 'west');
%     ['Ellipsoid fit: a = ', num2str(fit_ell), ', avg. sq. dist. = ', num2str(min(score_ell) / size(X, 2))], 'Location', 'west');
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, ['media/moss_chloronema_hyphoid_fit.png']);
close all;

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