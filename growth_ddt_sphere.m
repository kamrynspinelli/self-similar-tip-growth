step = 0.001;
gam = zeros(1,int32(1/step) + 1);
t = 0.00001;

global l0 r0
l0 = zeros(1,floor((pi/2 + 1) / step));
r0 = zeros(1,floor((pi/2 + 1) / step));
for i = 1:floor((pi/2 + 1) / step) % get the initial configuration
    l0(i) = step;
    r0(i) = init_outline(i*step);
end

for growth_inc = 1:500 % do a bunch of growth increments
    gam(1) = 1;
    for i = 1:floor((pi/2 + 1) / step)
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

z0 = zeros(1,floor((pi/2 + 1) / step));
z0(1) = 2;
a = 0;
for i=1:floor((pi/2 + 1) / step) - 1 % get the shape after growth
    z0(i+1) = z0(i) - sqrt(l0(i)^2 - (r0(i) - a)^2);
    a = r0(i);
end
plot(z0 - min(z0), r0); daspect([1 1 1]);
% plot(z0 - min(z0), r0); hold on; plot(z0 - min(z0), gam_trim); daspect([1 1 1]);

function f = init_outline(s)
    if s < pi/2
        f = sin(s);
    else
        f = 1;
    end
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

function df = init_doutline(s)
    if s < pi/2
        df = cos(s);
    else
        df = 0;
    end
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