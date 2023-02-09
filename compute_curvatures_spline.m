function [ks, ktheta] = compute_curvatures_spline(rv)
    N = size(rv,2);
    
    % Decide which point will separate the horizontally- and
    % vertically-oriented splines. This will happen when dr/dz = slope_cutoff.
    % slope_cutoff = -1/2;
    slope_cutoff = -1 + 0.001;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope >= slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    
    % left of cut point
    for index = 1:cut-1
        z = rv(1,index); % coordinates of the point
        r = rv(2,index);
        r_p = 3 * horizA(index) * z^2 + 2 * horizB(index) * z + horizC(index); % r'(z)
        r_pp = 6 * horizA(index) * z + 2 * horizB(index); % r''(z)
        sin_alpha = 1 / sqrt(1 + r_p^2); % local angle 
        ks(index) = r_pp / (1 + r_p^2)^(3/2);
        ktheta(index) = sin_alpha / r;
    end
    
    % cut point
    z = rv(1,cut); % coordinates of the point
    r = rv(2,cut);
    r_p = 3 * horizA(cut-1) * z^2 + 2 * horizB(cut-1) * z + horizC(cut-1); % r'(z)
    r_pp = 6 * horizA(cut-1) * z + 2 * horizB(cut-1); % r''(z)
    sin_alpha = 1 / sqrt(1 + r_p^2); % local angle 
    ks(cut) = r_pp / (1 + r_p^2)^(3/2);
    ktheta(cut) = sin_alpha / r;
    
    % right of cut point
    for index = cut+1:N-1
        z = rv(1,index); % coordinates of the point
        r = rv(2,index);
        z_p = 3 * vertA(index-cut) * (-r)^2 + 2 * vertB(index-cut) * (-r) + vertC(index-cut); % r'(z)
        z_pp = 6 * vertA(index-cut) * (-r) + 2 * vertB(index-cut); % r''(z)
        sin_alpha = z_p / sqrt(1 + z_p^2); % local angle 
        ks(index) = z_pp / (1 + z_p^2)^(3/2);
        ktheta(index) = sin_alpha / r;
    end
    
    % tip point
    z = rv(1,N); % coordinates of the point
    r = rv(2,N);
    z_p = 3 * vertA(N-cut) * (-r)^2 + 2 * vertB(N-cut) * (-r) + vertC(N-cut); % r'(z)
    z_pp = 6 * vertA(N-cut) * (-r) + 2 * vertB(N-cut); % r''(z)
    ks(N) = z_pp / (1 + z_p^2)^(3/2);
    z_p_lim = 3 * vertA(N-cut) * (-rv(2,N-1)/100)^2 + 2 * vertB(N-cut) * (-rv(2,N-1)/100) + vertC(N-cut); % r'(z)
    sin_alpha = z_p_lim / sqrt(1 + z_p_lim^2); % local angle 
    % ktheta(N) = ktheta(N-1); % assume continuity
    ktheta(N) = sin_alpha / (rv(2,N-1)/100);
    
    ks = -ks; % make both curvatures positive for a sphere
end