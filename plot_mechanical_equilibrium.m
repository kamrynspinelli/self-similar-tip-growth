sS = sigmaS(strainl, strainr, k, mu);
sTheta = sigmaTheta(strainl, strainr, k, mu);
plot(s, fliplr(ksspline(2:end) .* sS(2:end) + kthetaspline(2:end) .* sTheta(2:end)), ...
    s, fliplr(kthetaspline(2:end) .* sS(2:end)), 'LineWidth', 2.0);
xlim([0 max(s)]);
ylim([0.4 1.1]);
pbaspect([1 1 1]);
xlabel('s');
legend('\kappa_s \sigma_s + \kappa_\theta \sigma_\theta', ...
    '\kappa_\theta \sigma_s', 'Location', 'east');
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, plotname);