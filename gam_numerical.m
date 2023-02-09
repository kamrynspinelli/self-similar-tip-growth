function gam = gam_numerical(strainl, strainr, L0Spline, R0SplineLocal, N)
    s(1) = 0; % first entry is the tip
    s = [0 cumsum(fliplr(L0Spline'))];
    gam(1) = 1;
    % gam_s(1) = 1;
    % gam_theta(1) = 1;
    gam_s(1) = strainl(end) - 1;
    gam_theta(1) = strainl(end) - 1;
    for i = 1:size(L0Spline, 1)-1
        gam_theta(i+1) = (R0SplineLocal(N-i+1) - R0SplineLocal(N-i+2)) / L0Spline(N-i+1) / R0SplineLocal(N-i+1) ...
            * gam_s(1:i) * L0Spline(N:-1:N-i+1);
        gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
        gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
    end
end