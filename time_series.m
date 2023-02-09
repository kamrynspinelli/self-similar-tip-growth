hold on;
plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
for frameNum = 1:20:(1+20*10)
    p = reshape(f(frameNum,:), N+1, 2)';
    % plot(p(1,:), R0f(frameNum,:), 'LineWidth', 2.0); % plot the intrinsic R0 data against the z data for this frame
    set(gca,'ColorOrderIndex',1); % reset the coloring for the rest of the plots
    slope_cutoff = -1 + 0.001; % get the spline data
    % slope_cutoff = -0.3 + 0.001; % get the spline data
    cut = 2;
    slope = (p(2,3) - p(2,1)) / (p(1,3) - p(1,1));
    while slope > slope_cutoff
        cut = cut+1;
        slope = (p(2,cut+1) - p(2,cut-1)) / (p(1,cut+1) - p(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(p(1,1:cut), p(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-p(2,cut:end), p(1,cut:end), -1/slope, 0);
    for i = 1:size(horizA,2)
        spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
        plot(p(1,i):0.001:p(1,i+1), spl(p(1,i):0.001:p(1,i+1)), 'LineWidth', 2.0);
    end
    for i = 1:size(vertA,2)
        spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
        plot(spl(-p(2,cut+i-1):0.001:-p(2,cut+i)), -(-p(2,cut+i-1):0.001:-p(2,cut+i)), 'LineWidth', 2.0);
    end
    % xlim([0 7]); ylim([0 1.5]);
    xlim([0 rv(1,end)+0.5]); ylim([0 1.5]);
    % pause(0.2);
end
daspect([1 1 1]);
set(gca,'visible','off');
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% xlim([3 rv(1,end)+0.1])