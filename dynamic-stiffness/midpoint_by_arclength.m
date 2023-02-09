function mdpt = midpoint_by_arclength(df, x1, x2)
    % midpoint_by_arclength
    %   df: derivative of the function f:R->R
    %   x1, x2: real numbers with x1<x2
    %   returns: the x-value st. the arclength of f from x1 to x is
    %   equal to the arclength from x to x2
	tol = 0.001; % tolerance in arclength at the computed point
    % tol = 0.00001; % tolerance in arclength at the computed point
    ds = @(x) sqrt(1+(df(x)).^2); % the arclength element of f
%     segment_arclength = integral(ds, x1, x2, 'ArrayValued', true); % built-in numerical integration
	% segment_arclength = midpt(ds, x1, x2, 16); % midpoint rule
    if x1 > x2 % need to swap these to avoid a negative integral
        tmp = x1;
        x1 = x2;
        x2 = tmp;
    end
    segment_arclength = midpt(ds, x1, x2, 256); % midpoint rule
    % arclength from x1 to x2
    half_arclength = segment_arclength/2;
    mdpt = (x1+x2)/2;
%     while abs(integral(ds,x1,mdpt) - half_arclength) > tol % Newton's
%         % method, built-in numerical integration
%         mdpt = mdpt - (integral(ds, x1, mdpt,'ArrayValued', true) - half_arclength)/ds(mdpt);
%     end
%     while abs(midpt(ds, x1, mdpt, 16) - half_arclength) > tol % Newton's method, midpoint rule integration
%         mdpt = mdpt - (midpt(ds, x1, mdpt, 16) - half_arclength)/ds(mdpt);
%     end
    while abs(midpt(ds, x1, mdpt, 256) - half_arclength) > half_arclength * tol % Newton's method, midpoint rule integration
        mdpt = mdpt - (midpt(ds, x1, mdpt, 256) - half_arclength)/ds(mdpt);
    end
%     l = x1; % bounds for bisection method
%     r = x2;
%     mdpt = (l+r)/2;
%     while abs(midpt(ds, x1, mdpt, 256) - half_arclength) > tol % Bisection method
%         if (midpt(ds, x1, mdpt, 256) > half_arclength && x1 < x2) || ...
%                 (midpt(ds, x1, mdpt, 256) < half_arclength && x1 > x2)
%             r = mdpt;
%             mdpt = (l+r)/2;
%         else
%             l = mdpt;
%             mdpt = (l+r)/2;
%         end
%     end
end