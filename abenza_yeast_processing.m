clear all;
load('../cell-profiles/abenza-yeast/actual/cropreg347-01_03_R3D_PRJ-t.mat'); % load the data

N = 64; % the desired number of points for each outline

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

% do the rotation
angle = - atan(-1/a) + pi/2;
rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
for i = 1:size(Cx,2)
    rottmp = rot * [Cx(:,i) Cy(:,i)-b]';
    Cx_rot(:,i) = rottmp(1,:);
    Cy_rot(:,i) = rottmp(2,:);
end
hold on; for i = 1:size(Cx,1); plot(Cx_rot(:,i), Cy_rot(:,i)); end; daspect([1 1 1]);
exportgraphics(gcf, 'media/abenza-yeast-processing-after-rotation.png');
close all;

% make discretized outlines
% top half
kssumupper = zeros(1,N);
kthetasumupper = zeros(1,N);

hold on;
% for prof = 1:10
for prof = 1:size(Cx,1)
    X = Cx_rot(1:32,prof)'; % load the x and y coordinates for this profile. 1:32 selects the old end of the yeast
    Y = Cy_rot(1:32,prof)'; % they trace the profile from bottom left, to tip, to top left
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
    
    [ks ktheta] = compute_curvatures(rv);
    ksuppers{prof} = ks;
    kthetauppers{prof} = ktheta;
    
    kssumupper = kssumupper + ks;
    kthetasumupper = kthetasumupper + ktheta;
    
%     plot(0:N-1, fliplr(ks), 0:N-1, fliplr(ktheta));
end

% bottom half
kssumlower = zeros(1,N);
kthetasumlower = zeros(1,N);
% for prof = 1:10
for prof = 1:size(Cx,1)
    X = Cx_rot(1:32,prof)'; % load the x and y coordinates for this profile. 1:32 selects the old end of the yeast
    Y = Cy_rot(1:32,prof)'; % they trace the profile from bottom left, to tip, to top left
    Y = -Y; % use the lower profile instead of the top
    X = fliplr(X);
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
    rv(2,:) = -rv(2,:); % use the lower profile instead of the top
    rv(2,:) = rv(2,:) - min(rv(2,:)); % center at (0,0)
    rv(1,:) = rv(1,:) - min(rv(1,:));
    rv = rv / max(rv(2,:)); % normalize to have height 1
    rv_lowers{prof} = rv;
    
    [ks ktheta] = compute_curvatures(rv);
    kslowers{prof} = ks;
    kthetalowers{prof} = ktheta;
    kssumlower = kssumlower + ks;
    kthetasumlower = kthetasumlower + ktheta;
    plot(-N+1:0, ks, -N+1:0, ktheta);
end
close all;

% plot(-N:0, kssumlower / size(Cx,1), -N:0, kthetasumlower / size(Cx,1), 0:N, fliplr(kssumupper) / size(Cx,1), 0:N, fliplr(kthetasumupper) / size(Cx,1))
hold on;
% for i = 1:size(Cx,1)
% for i = [1 7 17 18 19 20 21 23 24 26] % 10 profiles where the curvatures match nicely at the tip
% %     close all; hold on;
% %     plot(-N:0, kslowers{i}, 0:N, fliplr(ksuppers{i}));
%     plot(-N:0, kthetalowers{i}, 0:N, fliplr(kthetauppers{i}));
% end
inds = [1 7 17 18 19];
% for i = inds % just use 5 of the profiles but treat the upper and lower parts as separate outlines
% %     close all; hold on;
%     plot(0:N-1, fliplr(kslowers{i}), 0:N-1, fliplr(ksuppers{i}));
% %     plot(0:N-1, fliplr(kthetalowers{i}), 0:N-1, fliplr(kthetauppers{i}));
% end
% exportgraphics(gcf, 'media/abenza-yeast-ks-aggregate.png');
% exportgraphics(gcf, 'media/abenza-yeast-ktheta-aggregate.png');
close all;
for i = 1:5
    rvs_export{2 * i - 1} = rvs{i};
    rvs_export{2* i} = rv_lowers{i};
end

rvs = rvs_export;
save('../cell-profiles/abenza-yeast/actual/processed.mat', 'rvs')