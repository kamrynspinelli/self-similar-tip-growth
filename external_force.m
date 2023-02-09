function force = external_force(rv)
% Computes the effect of external forces on the membrane.
%   Inputs:
%       rv: the 2xn matrix of positions of points on the membrane
%   Returns:
%       force: a 2xn matrix holding the external forces on each point in 
%           the x and y directions
    global Fext external_force_mode variance_factor triangle_half_width;
    force = zeros(size(rv)); % create a zero matrix to hold the forces
    if external_force_mode == "tip_only"  % for tip_only mode
        force(1,end) = -Fext; % apply the external force only to the tip
    elseif external_force_mode == "gaussian" % for gaussian mode
        % calculate variance of the y-values of each point from the mean 0 (as 
        % the membrane is symmetric about the x-axis), where each point is
        % given probability 1/(number of points)
        yvariance = (1/length(rv))*sum(rv(2,:).^2) * variance_factor; % the variance
        sumoverpoints = 0; % this will be the sum of the values of the PDF at the y-value of each point
        for point = 1:length(rv)
            sumoverpoints = sumoverpoints + normpdf(rv(2,point), 0, sqrt(yvariance)); % accumulate the PDF value
        end
        for point = 2:length(rv) % exclude the leftmost point to avoid self-intersections of the membrane
            force(1,point) = -(Fext * normpdf(rv(2,point), 0, sqrt(yvariance)) / sumoverpoints);
        end
    elseif external_force_mode == "triangular"
        % imagine a right triangle whose right angle is at the cell tip,
        % opening up and to the left
        h = triangle_half_width; % height of the triangle
        F = Fext; % total area of the triangle
        top_point = find(rv(2,:) < h, 1);
        num_points = length(rv) - top_point + 1;
        b = -F / (sum(rv(2,top_point:end)) / h - num_points);
        for point = 2:length(rv)
            force(1,point) = min(0, -b + b / h * rv(2,point));
        end
    end
end

