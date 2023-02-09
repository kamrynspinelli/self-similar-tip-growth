function [L0 R0] = spline_intrinsics(rv)
    %% Find the splines
    % Decide which point will separate the horizontally- and
    % vertically-oriented splines. This will happen when dr/dz = ~slope_cutoff.
    N = size(rv, 2);
    % slope_cutoff = -1/2;
	slope_cutoff = -1 + 0.001;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope >= slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    
    % ===== DEBUG =====
%     hold on;
%     plot(rv(1,:), rv(2,:), 'o');
%     for i = 1:size(horizA,2)
%         spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
%         plot(rv(1,i):0.001:rv(1,i+1), spl(rv(1,i):0.001:rv(1,i+1)));
%     end
%     for i = 1:size(vertA,2)
%         spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
%         plot(spl(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), -(-rv(2,cut+i-1):0.001:-rv(2,cut+i)));
%     end
%     daspect([1 1 1]);
    % ===== END DEBUG =====
    
    %% Compute the intrinsic lengths and radii
    L0 = zeros(N-1,1);
    bisector = zeros(N-1,1);
    R0 = zeros(N-1,1);
    for i = 1:cut-1
        L0(i) = midpt(@(z) sqrt(1 + (3 * horizA(i) * z^2 + 2 * horizB(i) * z + horizC(i))^2), rv(1,i), rv(1,i+1), 256);
        bisector(i) = midpoint_by_arclength(@(z) 3 * horizA(i) * z.^2 + 2 * horizB(i) * z + horizC(i), rv(1,i), rv(1,i+1));
%         bisector(i) = (rv(1,i) + rv(1,i+1)) / 2;
        R0(i) = horizA(i) * bisector(i)^3 + horizB(i) * bisector(i)^2 + horizC(i) * bisector(i) + horizD(i);
    end
    for i = cut:N-1
        L0(i) = midpt(@(r) sqrt(1 + (3 * vertA(i-cut+1) * r^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1))^2), -rv(2,i), -rv(2,i+1), 256);
        bisector(i) = midpoint_by_arclength(@(r) 3 * vertA(i-cut+1) * r.^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1), -rv(2,i), -rv(2,i+1));
%         bisector(i) = (-rv(2,i) + -rv(2,i+1)) / 2;
        R0(i) = -bisector(i);
    end
end