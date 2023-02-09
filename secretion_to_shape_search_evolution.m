function [tip_r, elong_vel, gam_mass, steady_radius, strain_shank] =  secretion_to_shape_search_evolution(mean, std, init_N, num_steps, stiffness)
% gam is a function handlle with gam(0) = 1, init_N is the initial number
% of patches in the system, and num_steps is the number of iterations to
% run the growth and solver in order to reach the desired tip shape
global Tol Rtol TolFun TolX Inc P

gam = @(x) 1/(std * sqrt(2*pi)) * (exp(-1/2 * ((x-mean) / std).^2) + exp(-1/2 * ((x+mean) / std).^2));

arclength_cap = 15;

% stiffness = 24; % stiffness parameter for the deformation
N = init_N; % number of patches for the simulation
% tstep = 0.1; % timestep for the growth % TESTING
tstep = 0.5; % timestep for the growth % FOR SIMULATIONS

% initial guess based on hyphoid
a = 6.8; % hyphoid parameter for the initial guess
f = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
fp = @(y) pi/a .* cot(pi*y) - pi^2*y/a.*csc(pi*y).^2;
% find the y such that f(y) = -5, this will be the base of the cell
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

% ===== DEBUG =====
% rv0 = 4*rv0;
% ===== END DEBUG =====

P = 1; % set the pressure
k = stiffness; % set up the material property
mu = stiffness;
Rtol = 1e-13; % computational parameters
Tol = 1e-3;
TolFun = 1e-16;
TolX = 1e-16;
Inc =  0.25; % fast movie speed
% Tol = 3e-2;
% Inc = 0.25;

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

base_discretization = L0SplineLocal(2); % remeshing data

tic;

[ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
s = [0 cumsum(fliplr(L0Spline'))]; % indexed with tip first

for i = 1:num_steps    
%     if distortion(L0SplineLocal) >= 1.3
    if i > 1
        [N, LUT, rv, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal] = ...
            remesh(N, LUT, rv, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, base_discretization);
        L0Spline = diff(L0SplineLocal)';
    end
%     end
    
    if max(s) > arclength_cap % truncate the rear of the profile
        [rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, N] ...
            = trunc(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, N, arclength_cap);
        L0Spline = diff(L0SplineLocal)';
    end
    
    [rv, strainl, strainr] = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
    
    s = [0 cumsum(fliplr(L0Spline'))]; % indexed with tip first
    
    gam_s = gam(s) .* fliplr(strainl - 1); % \gamma_s(s) = \gamma(s) * (\lambda_s(s) - 1); indexed with tip first
    gam_theta = gam(s) .* fliplr(strainr - 1); % similarly for \gamma_\theta

%     L0SplineNew = L0Spline + tstep * fliplr(fliplr(L0Spline') .* gam_s(1:end-1))';
%     L0SplineLocal = [0 cumsum(L0Spline)'];
%     L0Spline = L0SplineNew;
    L0SplineLocalTip = -fliplr(L0SplineLocal - max(L0SplineLocal));
    for pt = 1:N+1
        % indexed with tip first
        L0SplineLocalTipNew(pt) = L0SplineLocalTip(pt) + tstep * sum(gam_s(1:pt-1) .* fliplr(L0Spline(end-pt+2:end)'));
    end
    L0SplineLocal = -fliplr(L0SplineLocalTipNew - max(L0SplineLocalTipNew));
    clearvars L0SplineLocalTipNew; % clear the variable to avoid problems after truncation
    L0Spline = diff(L0SplineLocal)';
    R0SplineLocal = R0SplineLocal + tstep * fliplr(gam_theta) .* R0SplineLocal;
    clf;
    plot_shape(rv);
    
    if i == 14
        1;
    end
end

toc;

[ksspline kthetaspline] = compute_curvatures_spline(rv);

gam_num = gam_numerical(strainl, strainr, L0Spline, R0SplineLocal, N);
tip_r = 1 / ksspline(end); % tip radius of curvature 
elong_vel = elongation_velocity(s, gam, strainl, rv);
gam_mass = gamma_mass(s, gam, strainl, rv);
steady_radius = steady_state_radius(rv);
strain_shank = strain_apical_shank(strainl, rv);
save(['secretion_to_shape_search_evolution_N_', num2str(init_N), '_frames_', num2str(num_steps), '_mean_', num2str(mean), '_std_', num2str(std), '_stiffness_', num2str(stiffness), '_truncated.mat']);
end

function [rv, strainl, strainr] = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun * 4,'TolX',TolX * 4,'Algorithm','levenberg-marquardt');
%     options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt');
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

function dis = distortion(A)
    % measures how nonuniform the mesh A is, by measuring the ratio of the 
    % size difference between largest and smallest patches
    dis = max(diff(A)) / min(diff(A));
end

function [rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, N] ...
    = trunc(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, N, arclength_cap)
    % truncates the rear of the cell profile, removing points from the
    % shank
    L0SplineLocalTip = -fliplr(L0SplineLocal - max(L0SplineLocal));
    % find the index to cut off at in order to retain the desired arclength
    % just below arclength_cap
    ind = N + 3 - min(find(L0SplineLocalTip > arclength_cap));
    rv = rv(:,ind:end);
    rv(1,:) = rv(1,:) - rv(1,1); % translate the profile to meet the r-axis
    N = size(rv,2)-1;
    L0SplineLocal = L0SplineLocal(ind:end);
    L0SplineLocal = L0SplineLocal - min(L0SplineLocal);
    R0SplineLocal = R0SplineLocal(ind:end);
    KSplineLocal = KSplineLocal(ind:end);
    MUSplineLocal = MUSplineLocal(ind:end);
    LUT = diag(ones(N+1,1)) + diag(-ones(N,1), 1);
    LUT = LUT(1:end-1,:);
end

function plot_shape(rv)
    set(gca,'ColorOrderIndex',1); % reset the coloring for the rest of the plots
    slope_cutoff = -1 + 0.001; % get the spline data
    % slope_cutoff = -0.3 + 0.001; % get the spline data
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope > slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    hold on; % do the plot
    for i = 1:size(horizA,2)
        spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
        plot(rv(1,i):0.001:rv(1,i+1), spl(rv(1,i):0.001:rv(1,i+1)), 'LineWidth', 2.0);
    end
    for i = 1:size(vertA,2)
        spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
        plot(spl(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), -(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), 'LineWidth', 2.0);
    end
    daspect([1 1 1]);
    xlim([0 max(rv(1,:))+0.5]); ylim([0 1.3]);
end

function ev = elongation_velocity(s, gam, strainl, rv)
    % s is stress-free arclength (really s^0)
%     ev = strainl(2) * sum(diff(s) .* gam(s(1:end-1)) .* fliplr(strainl(2:end) - 1));
    ind = max(find(diff(abs(diff(rv(2,:)) ./ diff(rv(1,:)))) < 0));
    if isempty(ind) % if the first node is the shank
        ind = 1;
    end
    ev = strainl(ind) * sum(diff(s(1:ind+1)) .* gam(s(1:ind)) .* fliplr(strainl(end-ind+1:end) - 1));
end

function m = gamma_mass(s, gam, strainl, rv)
%     m = sum(gam(s(1:end-1)) .* diff(s) .* fliplr(strainl(2:end)));
    ind = max(find(diff(abs(diff(rv(2,:)) ./ diff(rv(1,:)))) < 0));
    if isempty(ind) % if the first node is the shank
        ind = 1;
    end
    m = sum(gam(s(1:ind)) .* diff(s(1:ind+1)) .* fliplr(strainl(end-ind+1:end)));
end

function r = steady_state_radius(rv)
    ind = max(find(diff(abs(diff(rv(2,:)) ./ diff(rv(1,:)))) < 0));
    if isempty(ind) % if the first node is the shank
        ind = 1;
    end
    r = rv(2,ind);
end

function s = strain_apical_shank(strainl, rv)
    ind = max(find(diff(abs(diff(rv(2,:)) ./ diff(rv(1,:)))) < 0));
    if isempty(ind) % if the first node is the shank
        ind = 1;
    end
    s = strainl(ind);
end