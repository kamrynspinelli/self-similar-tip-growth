step = 0.01;
gam = zeros(1,int32(1/step) + 1);
t = 0.001;
a = 3;
b = 1;
quarter_arclength = ellipse_arclength(a, b, 0);

global l0 r0
l0 = zeros(1,floor((quarter_arclength + 1) / step));
r0 = zeros(1,floor((quarter_arclength + 1) / step));
z0 = zeros(1,floor((quarter_arclength + 1) / step));
for i = 1:floor((quarter_arclength + 1) / step) % get the initial configuration
    [z0(i), r0(i)] = init_outline(a, b, i*step);
end
l0(1) = sqrt((z0(1) - (1 + a))^2 + (r0(1))^2);
for i = 2:floor((quarter_arclength + 1) / step)
    l0(i) = sqrt((z0(i) - z0(i-1))^2 + (r0(i) - r0(i-1))^2);
end
clear z0; % this name will be used again later, get rid of it

for growth_inc = 1:200 % do a bunch of growth increments
    gam(1) = 1;
    for i = 1:floor((quarter_arclength + 1) / step)
        gam(i+1) = doutline(i*step + step * sum(t * gam(1:i))) ...
            / outline(i*step + step * sum(t * gam(1:i))) * (step * sum(gam(1:i)));
        if abs(gam(i+1)) > 5
            1;
        end
    end
    % plot(step * (0:floor((pi/2 + 1) / step)), gam); % s-gamma plot
    % plot(1 - cos(step * (0:floor((pi/2 + 1) / step))), gam); % x-gamma plot
    gam_trim = gam(2:end);
    l0 = l0 .* gam_trim * t + l0;
    r0 = r0 .* gam_trim * t + r0;
end

z0 = zeros(1,floor((quarter_arclength + 1) / step));
z0(1) = 1 + a;
tmp = 0;
for i=1:floor((quarter_arclength + 1) / step) - 1 % get the shape after growth
    z0(i+1) = z0(i) - sqrt(l0(i)^2 - (r0(i) - tmp)^2);
    tmp = r0(i);
end
plot(z0 - min(z0), r0); daspect([1 1 1]);
% plot(z0 - min(z0), r0); hold on; plot(z0 - min(z0), gam_trim); daspect([1 1 1]);

function [x, y] = init_outline(a, b, s)
    if s < ellipse_arclength(a, b, 0)
        % bisection method to find corresponding x
        ll = 0;
        ul = a;
        midpt = (ll + ul) / 2;
        tol = 0.01;
        while abs(ellipse_arclength(a, b, midpt) - s) > tol
            if ellipse_arclength(a, b, midpt) - s < 0 % didn't go far enough along
                ul = midpt;
            else
                ll = midpt;
            end
            midpt = (ll + ul) / 2;
        end
        x = midpt;
        y = b / a * sqrt(a^2 - x^2);
        x = x + 1; % account for the horizontal shift in initial profile
    else
        x = 1 - s + ellipse_arclength(a, b, 0);
        y = 1;
    end
end

function s = ellipse_arclength(a, b, x)
    s = midpt(@(t) sqrt(1 + b^2 / a^2 * t^2 / (a^2 - t^2)), x, a, 512);
end

function f = outline(s)
    global l0 r0
    i = 1;
    l_sum = l0(i);
    while l_sum < s && i < size(l0, 2)-1
        i = i+1;
        l_sum = l_sum + l0(i);
    end
    f = r0(i);
end

function df = doutline(s)
    global l0 r0
    i = 1;
    l_sum = l0(i);
    while l_sum < s && i < size(l0, 2)-1
        i = i+1;
        l_sum = l_sum + l0(i);
    end
    df = (r0(i+1) - r0(i)) / l0(i);
end