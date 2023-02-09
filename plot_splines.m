function [] = plot_splines(rv)
    %% Find the splines
    % Decide which point will separate the horizontally- and
    % vertically-oriented splines. This will happen when dr/dz = ~slope_cutoff.
    N = size(rv, 2);
    slope_cutoff = -1;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope > slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    
    %% Plot the splines
    hold on;
%     plot(rv(1,:), rv(2,:), 'o');
    for i = 1:size(horizA,2)
        spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
        plot(rv(1,i):0.001:rv(1,i+1), spl(rv(1,i):0.001:rv(1,i+1)), 'LineWidth', 2.0);
    end
    for i = 1:size(vertA,2)
        spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
        plot(spl(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), -(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), 'LineWidth', 2.0);
    end
    daspect([1 1 1]);
    hold off;
end