load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
indices = size(Px, 2) - [0 5 10 17 20 26 30 34 40 45];
f = figure;
hold on;
c = lines(8); % some colors to use
for profile_num = 1:5
    profile = indices(profile_num);
    X = Px(1:max(find(Px(:,profile))),profile)'; % load the x and y coordinates
    Y = Py(1:max(find(Px(:,profile))),profile)'; % they trace the profile from bottom left, to tip, to top left
%     Y = -Py(1:max(find(Px(:,profile))),profile)'; % use the lower profile instead of the top
    Y = Y - min(Y) - (max(Y) - min(Y))/2;
    X = X / max(abs(Y)) - 2.25;
    Y = Y / max(abs(Y));
    Indfirst = min(find(atan(diff(Y) ./ diff(X)) >= 0)); % indices of the boundary points - exclude points where the outline curves in at the shank
    Indlast = max(find(atan(diff(Y) ./ diff(X)) <= 0)) + 1;
    X = X(Indfirst:Indlast);
    Y = Y(Indfirst:Indlast);
%     Y = Y - min(Y) + 1;
    [AX BX CX DX] = find_splines(1:size(X,2), X, 1, -1); % interpolate the x and y coordinates with cubic splines
    [AY BY CY DY] = find_splines(1:size(X,2), Y, 0, 0);
    for i = 1:size(X,2)-1
        splX = @(x) AX(i) * x .^3 + BX(i) * x .^2 + CX(i) * x + DX(i);
        splY = @(x) AY(i) * x .^3 + BY(i) * x .^2 + CY(i) * x + DY(i);
        p = plot(splX(i:0.001:i+1), splY(i:0.001:i+1), 'LineWidth', 2.0);
        p.Color = c(profile_num,:);
    end
end
daspect([1 1 1]);
pbaspect([1 1 1]);
xlim([0 3.2]);
% set(gca,'visible','off');
% set(gca,'xticklabel',[],'yticklabel',[])
xlabel('z');
ylabel('r');
ax = gca;
set(ax, 'fontsize', 14);
% exportgraphics(gcf, 'media/dumais-root-hair-profiles.png');
f.PaperPosition = [1.333 3.31 5 5];
filename = 'media/dumais-root-hair-profiles.png';
print(filename, '-dpng', '-r100');
close all;