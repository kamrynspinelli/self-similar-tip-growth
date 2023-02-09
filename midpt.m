function int = midpt(f, xl, xr, div)
    % midpt
    %   Performs a numerical integration using the midpoint rule. Accepts a
    %   function handle f which is the function to be integrated, scalars
    %   xl and xr which are the endpoints of the interval of integration,
    %   and an integer div which is the number of subintervals the interval
    %   should be partitioned into for the integration. f must be scalar-
    %   or (column) vector-valued.
    step = (xr - xl) / div; % the length of the subintervals
    %part = linspace(xl, xr, div + 1); % the marker points of the partition
    mids = linspace(xl + step / 2, xr - step / 2, div); % the midpoints of each subinterval
    fvals(:,div) = f(mids(end)); % preallocate the values matrix with the last entry
    for i = 1:div-1 % and compute f at all the other midpoints
        fvals(:,i) = f(mids(i));
    end
    % add up all the values of f (summing over rows in the matrix) times the length of the subintervals
    int = step * sum(fvals, 2);
end