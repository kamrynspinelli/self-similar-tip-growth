function int = trapz_helper(f, xl, xr, div)
    step = (xr - xl) / div; % the length of the subintervals
    part = linspace(xl, xr, div + 1); % the marker points of the partition
    fvals(:,div+1) = f(part(end)); % preallocate the values matrix with the last entry
    for i = 1:div % and compute f at all the other midpoints
        fvals(:,i) = f(part(i));
    end
     % add up all the values of f (summing over rows in the matrix) times the length of the subintervals
    int = trapz(part, fvals, 2);
end