function int = simps(f, xl, xr, div)
    % simps
    %   Performs a numerical integration using Simpson's rule. Accepts a
    %   function handle f which is the function to be integrated, scalars
    %   xl and xr which are the endpoints of the interval of integration,
    %   and an integer div which is the number of subintervals the interval
    %   should be partitioned into for the integration. f must be scalar-
    %   or (column) vector-valued.
    step = (xr - xl) / div; % the length of the subintervals
    part = linspace(xl, xr, div + 1); % the marker points of the partition
    fvals(:,div+1) = f(part(end)); % preallocate the values matrix with the last entry
    for i = 1:div % and compute f at all the other midpoints
        fvals(:,i) = f(part(i));
    end
    int = fvals(:,1) + fvals(:,div+1);
    for i = 2:div
        int = int + (2 + 2*mod(i+1,2)) * fvals(:,i);
    end
    % add up all the values of f (summing over rows in the matrix) times the length of the subintervals
    int = step / 3 * int;
end