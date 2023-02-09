hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
ax.YColor = 'black';
xlim([0 3.2]);
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma
p4 = plot(rv(1,2:end), velocity, '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
p4.Color = c(5,:); % plot velocity
ax.YColor = 'black';
ylabel("secretion rate, velocity");
ax = gca;
ax.XAxis.Visible = 'off'; % turn off this x-axis
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
ax.YColor = 'black';
xlim([0 3.2]);
line([3.2 3.3], [0 0])
yyaxis right;
p5 = plot(rv(1,2:end), strainl, '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain \lambda_s');
p5.Color = c(3,:); % plot meridional strain
p6 = plot(rv(1,2:end), strainr, '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain \lambda_\theta');
p6.Color = c(4,:); % plot circumferential strain
ymin = min([strainl strainr]); % adjust the y bounds for the strain plot
ymax = max([strainl strainr]);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
ax.YColor = 'black';
ylabel("strains");
ax = gca;
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top