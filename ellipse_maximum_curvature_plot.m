% plots the arclength coordinate of maximum curvature on an ellipsoid as a
% function of the aspect ratio a/b, where the ellipse is given by 
% b/a * sqrt(a^2 - x^2)
clear all;
integrand = @(a,t) sqrt(1 + (-t/a .* (a^2 - t.^2).^(-1/2)).^2);
a = 0.1:0.02:1;
for i = 1:size(a,2)
    s(i) = midpt(@(t) integrand(a(i), t), 0, a(i), 65536);
end
a = [a 1:0.02:2];
s = [s zeros(1,size(1:0.02:2, 2))];
plot(a, s);
plot(a, s, 'LineWidth', 2.0);
xlabel('a/R_N'); ylabel('s_k/R_N');
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, 'media/ellipsoid-location-of-max-curvature.png');
close all;

% example curvature profiles
s = 0:0.001:2.5;
xrng = 0:0.001:1.2;
arng = [1.1 1 0.9 0.7 0.5];
hold on;
for ai = 1:size(arng,2)
    clear s curv_rng;
    f = @(x) 1/arng(ai) * sqrt(arng(ai)^2 - x.^2);
    fp = @(x) -x ./ arng(ai) .* (arng(ai)^2 - x.^2).^(-1/2);
    fpp = @(x) -arng(ai) * (arng(ai)^2 - x.^2).^(-3/2);
    curv = @(x) -fpp(x) ./ (1 + fp(x).^2).^(3/2);
    for i = 1:size(xrng,2)-1
        curv_rng(i) = curv(xrng(i));
        s(i) = midpt(@(x) sqrt(1 + fp(x).^2), xrng(i), arng(ai), 2048);
    end
    curv_rng = curv_rng(not(imag(s)));
    s = s(not(imag(s)));
    s(end) = 0;
    curv_rng(end) = arng(ai);
    plot(s, curv_rng, 'LineWidth', 2.0);
    ax = gca;
    set(gca,'ColorOrderIndex',ai+1); % reset the coloring for the rest of the plots
    if arng(ai) ~= 1
        plot(s(find(curv_rng == max(curv_rng))), curv_rng(find(curv_rng == max(curv_rng))), 'xk', 'LineWidth', 2.0); 
    end
end
xlabel('s');
ylabel('\kappa_s');
xlim tight;
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, 'media/ellipsoid-sample-curvature-profiles.png');
close all;