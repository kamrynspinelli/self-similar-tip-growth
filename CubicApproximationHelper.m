classdef CubicApproximationHelper
    % CubicApproximationHelper
    % Contains methods to help with the cubic approximation.
    methods (Static = true)
        function coeff = cubic_coefficients(points)
            % cubic_coefficients
            %   points: a 2x4 matrix whose columns are points of the form
            %   (x,y)
            %   returns: a vector [a b c d] such that ax^3+bx^2+cx+d is the
            %   unique cubic passing through all four points
            y = (points(2,1:4))'; % column vector of y-values
            m = vander(points(1,1:4)); % vandermonde matrix of the x-values
            coeff = (m\y)'; % the coefficients of the cubic
        end
        
        function coeff = all_cubic_coefficients(points)
            % all_cubic_coefficients
            %   points: a 2xn matrix, n>=3, whose columns are points of the
            %   form (x,y)
            %   returns: a 4x(n-2) matrix whose columns have the form 
            %   [a b c d]'; the ith column is the coefficients of the cubic
            %   approximant for the ith patch described by the given
            %   sequence of points
            n = length(points);
            ext_points = [points(:,2) points(:,:)]; % the extended matrix 
            % of points with the second point mirrored over the y-axis in 
            % the first entry
            ext_points(1,1) = -ext_points(1,1); % negate that x-value
            for i=1:n-2
                slices(i,:,:) = ext_points(:,i:i+3); % the points needed to
                % compute the cubic coefficients for the ith patch
                coeff(:,i) = CubicApproximationHelper.cubic_coefficients(reshape(slices(i,:,:), 2, 4));
                % the coefficients for the ith patch
            end            
        end
        
        function cubic = cubic_approximant(points)
            % cubic_approximant
            %   points: a 2x4 matrix whose columns are points of the form
            %   (x,y)
            %   returns: a function handle for the unique cubic passing
            %   through all four points
            coeff = CubicApproximationHelper.cubic_coefficients(points);
            % coefficients of the cubic
            cubic = @(x) coeff(1)*x.^3 + coeff(2)*x.^2 + coeff(3)*x ...
                + coeff(4); % the cubic with those coefficients
        end
        
        function mdpt = midpoint_by_arclength(df, x1, x2)
            % midpoint_by_arclength
            %   df: derivative of the function f:R->R
            %   x1, x2: real numbers with x1<x2
            %   returns: the x-value st. the arclength of f from x1 to x is
            %   equal to the arclength from x to x2
            tol = 0.001; % tolerance in arclength at the computed point
            ds = @(x) sqrt(1+(df(x)).^2); % the arclength element of f
            segment_arclength = integral(ds, x1, x2, 'ArrayValued', true);
            % arclength from x1 to x2
            half_arclength = segment_arclength/2;
            mdpt = (x1+x2)/2;
            while abs(integral(ds,x1,mdpt) - half_arclength) > tol % Newton's method
                mdpt = mdpt - (integral(ds, x1, mdpt,'ArrayValued', true) - half_arclength)/ds(mdpt);
            end
        end
    end
end

