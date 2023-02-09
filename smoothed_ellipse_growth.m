clear x s gam s_inc y;

a = 3;
b = 1;
xstep = 0.0001;
x = a+1:-xstep:0;
gam = zeros(1, floor((a+1) / xstep) + 1);
s = zeros(1, floor((a+1) / xstep) + 1);
s_inc = zeros(1, floor((a+1) / xstep) + 1);
y = zeros(1, floor((a+1) / xstep) + 1);
gam(1) = 1;
s(1) = 0;
y(1) = 0;
for i = 2:size(x,2)
    [s(i), y(i)] = outline(a, b, x(i));
    s_inc(i-1) = s(i) - s(i-1);
    gam(i) = ((y(i) - y(i-1)) / s_inc(i-1)) / y(i) * (s_inc(1:i-1) * gam(1:i-1)');
end
plot(s, gam, 'LineWidth', 2.0);
title(['Isotropic self-similar growth profile: smoothed C^2 ellipse, a = ', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([0 1]);
xlabel("s");
ylabel("\gamma(s)");
exportgraphics(gcf, ['media/isotropic_growth_c2_ellipse_a_', num2str(a), '_b_', num2str(b), '.png']);

function [s, y] = outline(a, b, x)
    if x >= 1
        s = arclength(a, b, x-1);
        y = 1 + (b / a * sqrt(a^2 - (x-1)^2) -1) * (6*(x-1)^5-15*(x-1)^4+10*(x-1)^3);
    else
        s = (1 - x) + arclength(a, b, 0);
        y = 1;
    end
end

function s = arclength(a, b, x)
    s = midpt(@(t) sqrt(1 + ((t^2 * (30 * a^2 * b * (-1 + t)^2 - 30 * a * (-1 + t)^2 * sqrt((a - t) * (a + t)) + b * t^2 * (-40 + 75 * t - 36 * t^2)))/(a * sqrt((a - t) * (a + t))))^2), x, a, 512);
end