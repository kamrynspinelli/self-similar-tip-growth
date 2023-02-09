clear x s gam s_inc y;

xstep = 0.0001;
x = 2:-xstep:0;
gam = zeros(1, floor((2) / xstep) + 1);
s = zeros(1, floor((2) / xstep) + 1);
s_inc = zeros(1, floor((2) / xstep) + 1);
y = zeros(1, floor((2) / xstep) + 1);
gam(1) = 1;
s(1) = 0;
y(1) = 0;
for i = 2:size(x,2)
    [s(i), y(i)] = outline(x(i));
    s_inc(i-1) = s(i) - s(i-1);
    gam(i) = ((y(i) - y(i-1)) / s_inc(i-1)) / y(i) * (s_inc(1:i-1) * gam(1:i-1)');
end
plot(s, gam, 'LineWidth', 2.0);
title(['Isotropic self-similar growth profile: smoothed C^2 sphere']);
xlim([0 max(s)]);
ylim([0 1]);
xlabel("s");
ylabel("\gamma(s)");
exportgraphics(gcf, ['media/isotropic_growth_c2_sphere.png']);

function [s, y] = outline(x)
    if x >= 1
        s = arclength(x-1);
        y = 1 + (sqrt(1-(x-1)^2)-1) * (6*(x-1)^5-15*(x-1)^4+10*(x-1)^3);
    else
        s = (1 - x) + arclength(0);
        y = 1;
    end
end

function s = arclength(x)
    s = midpt(@(t) sqrt(1 + (-((t * (10 * t^3 - 15 * t^4 + 6 * t^5))/sqrt(1 - t^2)) + (30 * t^2 - 60 * t^3 + 30 * t^4) * (-1 + sqrt(1 - t^2)))^2), x, 1, 512);
end