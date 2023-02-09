function secretion_to_shape_search(gam) % gam is a function handlle with gam(0) = 1
L2Tol = 0.01; % L^2 norm tolerance for gamma
stepTol = 1;
Del = 0.001; % distance increment for point perturbations
stiffness = 24; % stiffness parameter for the deformation
N = 48; % number of patches for the simulation

% initial guess based on secretion - fails
% s_max = 5; % assume the secretion is negligible past s = 5
% total_gam = midpt(gam, 0, s_max, 1024); % find the total amount of secretion
% rv0(1,1) = 0; % set up the initial outline, starting with the tip point 
% rv0(2,1) = 0;
% for i = 1:N % fill in the rest of the points
%     rv0(2,i+1) = 1/total_gam * midpt(gam, 0, s_max * i / N, 1024);
%     rv0(1,i+1) = rv0(1,i) - (s_max / N) * sqrt(1 - (1/total_gam * gam(s_max / N))^2);
% end
% rv0 = fliplr(rv0); % flip and shift to match model assumptions
% rv0(1,:) = rv0(1,:) - min(rv0(1,:));
% % we end up with a sharp point here - is it OK for the initial guess?

% initial guess based on hyphoid
a = 6.8; % hyphoid parameter for the initial guess
f = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
fp = @(y) pi/a .* cot(pi*y) - pi^2*y/a.*csc(pi*y).^2;
% find the y such that f(y) = -3, this will be the base of the cell
yl = 0; % bisection method
yr = 0.5;
while f(yr) > -5
    yr = (1 + yr) / 2;
end
yr_init = yr;
mid = (yl+yr)/2;
fmid = f(mid);
while abs(fmid - -5) > 0.001
    if fmid < -5
        yr = mid;
    else
        yl = mid;
    end
    mid = (yl+yr)/2;
    fmid = f(mid);
end
total_arclength = midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048);
% uniform point distribution
for i = 0:N
    target_arclength = total_arclength * i / N;
    yl = 0;
    yr = yr_init;
    mid = (yl+yr)/2;
    fmid = f(mid);
    while abs(midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) - target_arclength) > total_arclength / N * 0.001
        if midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) > target_arclength
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f(mid);
    end
    rv0(1,N+1-i) = f(mid);
    rv0(2,N+1-i) = mid;
end
rv0(1,N+1) = 1/a;
rv0(2,N+1) = 0;
rv0(1,:) = rv0(1,:) - min(rv0(1,:));

[s_turgid, gam_turgid, rv] = secretion(rv0, stiffness, N); % error associated with initial guess
dist = d(gam, s_turgid, gam_turgid);

for i = 1:N+1 % for each point, starting with the tip
    stepDist = Inf; % just to get things started
    stepDecreasing = 1;
    while (dist > L2Tol) && (stepDist > stepTol) && stepDecreasing % as long as the error is large and the step size hasn't decreased too much:
        % compute the gradient associated with perturbing the i-th marker point
        disp('Optimizing...');
        if i ~= N+1
            [s_turgid, gam_turgid, rv] = secretion(rv0 + Del * [zeros(2,N-i+1) [1; 0] zeros(2,i-1)], stiffness, N); % perturb the x coordinate a bit
            dx = (d(gam, s_turgid, gam_turgid) - dist) / Del;
        else
            dx = 0;
        end
        if i ~= 1
            [s_turgid, gam_turgid, rv] = secretion(rv0 + Del * [zeros(2,N-i+1) [0; 1] zeros(2,i-1)], stiffness, N); % perturb the y coordinate a bit
            dy = (d(gam, s_turgid, gam_turgid) - dist) / Del;
        else
            dy = 0;
        end
%         if sqrt(dx^2 + dy^2) > stepDist % if the step distance is increasing again, move on to the next point
%             stepDecreasing = 0;
%             disp(['Point: ', num2str(i), '; Dist: ', num2str(dist), '; stepDist: ', num2str(stepDist), ' (step length started increasing)']);
%         else
            stepDist = sqrt(dx^2 + dy^2);
            grad = [dx; dy];
            rv0 = rv0 + [zeros(2,N-i+1) -Del*grad zeros(2,i-1)]; % perturb the i-th marker point according to the gradient
            [s_turgid, gam_turgid, rv] = secretion(rv0, stiffness, N); % update the L2-distance
            dist = d(gam, s_turgid, gam_turgid);
            disp(['Point: ', num2str(i), '; Dist: ', num2str(dist), '; stepDist: ', num2str(stepDist)]);
%         end
    end
    
end

1;

end

function [s_turgid, gam_turgid, rv] = secretion(rv0, stiffness, N)
    global Tol Rtol TolFun TolX Inc P
    P = 1; % set the pressure
    k = stiffness; % set up the material property
    mu = stiffness;
    Rtol = 1e-13; % computational parameters
    Tol = 1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    Inc =  0.25; % fast movie speed
    
    rv = rv0;
    
    % make the LUT and ext_verts - I should really rework the code so these
    % aren't necessary anymore
    LUT = diag(ones(N+1,1)) + diag(-ones(N,1), 1);
    LUT = LUT(1:end-1,:);
    ext_verts = [1 N+1];
    
    initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
    for i=N:-1:1  % set the unturgid data
       index = find(LUT(i,:));
       L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs for ODE solver
       R0Parabolic(i) = initialArcs(i+1).vert(2);
    end    
    [L0Spline R0Spline] = spline_intrinsics(rv); % splines for nonlinear solver
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal = rv(2,:);
    K = k*ones(N,1);
    MU = mu*ones(N,1);  
    KSplineLocal = k * ones(N+1, 1);
    MUSplineLocal = mu * ones(N+1, 1);
    
    [rv, strainl, strainr] = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
%     [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % run the simulation
    
%     strainl = strainFramesl1(end,:);
%     strainr = strainFramesr1(end,:);
    s = zeros(1, size(L0Spline, 2) + 1);
    gam_turgid = zeros(1, size(L0Spline, 2) + 1);
    gam_s = zeros(1, size(L0Spline, 2) + 1);
    gam_theta = zeros(1, size(L0Spline, 2) + 1);
    s0(1) = 0; % first entry is the tip
    gam_turgid(1) = 1;
    gam_s(1) = strainl(end) - 1;
    gam_theta(1) = strainl(end) - 1;
    for i = 1:size(L0Spline, 1)-1
        s0(i+1) = s0(i) + L0Spline(N-i+1);
        gam_theta(i+1) = (rv0(2,N-i+1) - rv0(2,N-i+2)) / L0Spline(N-i+1) / rv0(2,N-i+1) ...
            * gam_s(1:i) * L0Spline(N:-1:N-i+1);
        gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
        gam_turgid(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
    end
    s_turgid = [0 cumsum(fliplr(L0Spline(2:end)' .* strainl(3:end)))]; % post-deformation arclength
end

function [rv, strainl, strainr] = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun * 4,'TolX',TolX * 4,'Algorithm','levenberg-marquardt','Display','off');
    X = fsolve(@solver_cspline_local_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,0); % do the deformation
    rv = reshape(X',N,2)'; 
    
    % compute the strains
    % Decide which point will separate the horizontally- and
    % vertically-oriented splines. This will happen when dr/dz = slope_cutoff.
    % slope_cutoff = -2;
    N = size(rv,2);
    slope_cutoff = -1 + 0.001;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope >= slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    arclength = zeros(N,1);
    arclength(1) = 0;
    for i = 1:cut-1
        arclength(i+1) = midpt(@(z) sqrt(1 + (3 * horizA(i) * z^2 + 2 * horizB(i) * z + horizC(i))^2), rv(1,i), rv(1,i+1), 256);
    end
    for i = cut:N-1
        arclength(i+1) = midpt(@(r) sqrt(1 + (3 * vertA(i-cut+1) * r^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1))^2), -rv(2,i), -rv(2,i+1), 256);
    end

    S = cumsum(arclength)'; % pairs of corresponding s and s^0 coordinates
    % S0 = cumsum(L0);
    S0 = L0;
    % generate splines to interpolate the increasing function s(s^0)
    % [SA SB SC SD] = find_splines(S0, S, S(2) / S0(2), 1 ,1);
    [SA SB SC SD] = find_splines(S0, S, S(2) / S0(2), (S(end) - S(end-1)) / (S0(end) - S0(end-1)));
    [RA RB RC RD] = find_splines(S, rv(2,:), 0, 1);
    [R0A R0B R0C R0D] = find_splines(S0, R0, 0, 1);
    for i = 1:N-1
        strainl(i) = 3 * SA(i) * S0(i)^2 + 2 * SB(i) * S0(i) + SC(i); % ds/ds^0
    end
    strainl(N) = 3 * SA(end) * S0(end)^2 + 2 * SB(end) * S0(end) + SC(end);
    strainr(1:N-1) = rv(2,1:end-1) ./ R0(1:end-1);
    strainr(N) = strainl(N); % by L'Hopital's rule
end

function dist = d(gam, s_turgid, gam_turgid)
    % gam is a function handle, gam_turgid is a discretized function
    % corresponding to the entries of s_turgid. Returns the L^2 norm of
    % (gam - gam_turgid).
    dist = sum((gam_turgid(1:end-1) - gam(s_turgid(1:end-1))).^2 .* diff(s_turgid));
end