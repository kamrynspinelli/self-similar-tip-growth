function [NNew, LUTNew, rvNew, L0SplineLocalNew, R0SplineLocalNew, KSplineLocalNew, MUSplineLocalNew] = remesh(N, LUT, rv, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, base_discretization)
    % The remeshing will take place if the increments of L0SplineLocal
    % become too uneven.
    NNew = ceil(max(L0SplineLocal) / base_discretization);
    LUTNew = zeros(NNew, NNew+1);
    for i = 1:NNew
        LUTNew(i,i) = 1;
        LUTNew(i,i+1) = -1;
    end
    
    % remesh L0 to uniform spacing, preserving the number of points
%     L0SplineLocalNew = linspace(0, max(L0SplineLocal), NNew+1); % uniform mesh in L0
    L0SplineLocalNew = fliplr([max(L0SplineLocal):-base_discretization:0 0]);
    
    % remesh R0
    R0SplineLocalNew = zeros(size(R0SplineLocal));
%     [R0A R0B R0C R0D] = find_splines(L0SplineLocal, R0SplineLocal, 0, -1);
    [R0A R0B R0C R0D] = find_splines(L0SplineLocal, R0SplineLocal, ...
        (R0SplineLocal(2)-R0SplineLocal(1))/(L0SplineLocal(2)-L0SplineLocal(1)), -1);
    for i = 1:size(L0SplineLocalNew, 2)
        ind = min(max(find(L0SplineLocal <= L0SplineLocalNew(i))), size(R0A, 2));
        R0SplineLocalNew(i) = R0A(ind) * L0SplineLocalNew(i)^3 ...
            + R0B(ind) * L0SplineLocalNew(i)^2 ...
            + R0C(ind) * L0SplineLocalNew(i) ...
            + R0D(ind);
    end
    R0SplineLocalNew(end) = 0; % don't let any numerical error in here!
    
    % remesh K
    KSplineLocalNew = zeros(size(KSplineLocal));
    [KA KB KC KD] = find_splines(L0SplineLocal, KSplineLocal, ...
        (KSplineLocal(2)-KSplineLocal(1)) / (L0SplineLocal(2)-L0SplineLocal(1)), ...
        (KSplineLocal(end)-KSplineLocal(end-1)) / (L0SplineLocal(end)-L0SplineLocal(end-1)));
    for i = 1:size(L0SplineLocalNew, 2)
        ind = min(max(find(L0SplineLocal <= L0SplineLocalNew(i))), size(KA, 2));
        KSplineLocalNew(i) = KA(ind) * L0SplineLocalNew(i)^3 ...
            + KB(ind) * L0SplineLocalNew(i)^2 ...
            + KC(ind) * L0SplineLocalNew(i) ...
            + KD(ind);
    end
    
    % remesh mu
    MUSplineLocalNew = zeros(size(MUSplineLocal));
    [MUA MUB MUC MUD] = find_splines(L0SplineLocal, MUSplineLocal, ...
        (MUSplineLocal(2)-MUSplineLocal(1)) / (L0SplineLocal(2)-L0SplineLocal(1)), ...
        (MUSplineLocal(end)-MUSplineLocal(end-1)) / (L0SplineLocal(end)-L0SplineLocal(end-1)));
    for i = 1:size(L0SplineLocalNew, 2)
        ind = min(max(find(L0SplineLocal <= L0SplineLocalNew(i))), size(KA, 2));
        MUSplineLocalNew(i) = MUA(ind) * L0SplineLocalNew(i)^3 ...
            + MUB(ind) * L0SplineLocalNew(i)^2 ...
            + MUC(ind) * L0SplineLocalNew(i) ...
            + MUD(ind);
    end
    
%     hold on; plot(L0SplineLocal, R0SplineLocal, '-', L0SplineLocalNew, R0SplineLocalNew, '-');
    
    % form the cubic spline system for the marker points
    slope_cutoff = -1 + 0.001;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope >= slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    
%     % get the correspondence between arclength and x/y-coordinates
%     arclengths = zeros(size(rv,2)-1,1);
%     for i = 1:cut-1
%         arclengths(i) = midpt(@(z) sqrt(1 + (3 * horizA(i) * z^2 + 2 * horizB(i) * z + horizC(i))^2), rv(1,i), rv(1,i+1), 256);
%     end
%     for i = cut:size(rv,2)-1
%         arclengths(i) = midpt(@(r) sqrt(1 + (3 * vertA(i-cut+1) * r^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1))^2), -rv(2,i), -rv(2,i+1), 256);
%     end
%     S = [0 cumsum(arclengths)];
%     [SL0A SL0B SL0C SL0C] = find_splines(L0SplineLocal, S, ...
%         (S(2) - S(1)) / (L0SplineLocal(2) - L0SplineLocal(1)), ...
%         (S(end) - S(end-1)) / (L0SplineLocal(end) - L0SplineLocal(end-1)));

    % get the correspondence between L0 and z/r-coordinates
    [ZL0A ZL0B ZL0C ZL0D] = find_splines(L0SplineLocal, rv(1,:), ...
        (rv(1,2) - rv(1,1)) / (L0SplineLocal(2) - L0SplineLocal(1)), 0);
    [RL0A RL0B RL0C RL0D] = find_splines(L0SplineLocal, rv(2,:), ...
        0, (rv(2,end) - rv(2,end-1)) / (L0SplineLocal(end) - L0SplineLocal(end-1)));
%     [ZL0A ZL0B ZL0C ZL0D] = find_splines(L0SplineLocal, rv(1,:), ...
%         1, 0);
%     [RL0A RL0B RL0C RL0D] = find_splines(L0SplineLocal, rv(2,:), ...
%         0, -1);
    
    % remesh rv
    rvNew = zeros(size(rv));
    for i = 1:size(L0SplineLocalNew, 2)
        ind = min(max(find(L0SplineLocal <= L0SplineLocalNew(i))), size(ZL0A, 2));
        rvNew(1,i) = ZL0A(ind) * L0SplineLocalNew(i)^3 ...
            + ZL0B(ind) * L0SplineLocalNew(i)^2 ...
            + ZL0C(ind) * L0SplineLocalNew(i) ...
            + ZL0D(ind);
        rvNew(2,i) = RL0A(ind) * L0SplineLocalNew(i)^3 ...
            + RL0B(ind) * L0SplineLocalNew(i)^2 ...
            + RL0C(ind) * L0SplineLocalNew(i) ...
            + RL0D(ind);
    end
end

