clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

name = 'hyphoid';
a = 6.8;
stiffness = 24;

tstep = 0.005;
num_frames = 2000;

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 128; %initial number of element

% Bulk Modulus
k = stiffness;
% Shear Modulus
mu = stiffness;
% Uniform material property
K = k*ones(N,1);
MU = mu*ones(N,1);

%computational parameters
Rtol = 1e-13;
Tol =  1e-3;
TolFun = 1e-16;
TolX = 1e-16;
% Inc =  1; % fast movie speed
% Inc = 0.25;
Inc = 0.05; % slow movie speed


%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated
L0Linear = ones(N,1);
R0Linear = ones(N,1);
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);

[ tmp, adj ] = generate_sphere(R,N); % get the adjacency data for N patches

f_shape = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
fp = @(y) pi/a .* cot(pi*y) - pi^2*y/a.*csc(pi*y).^2;

% find the y such that f(y) = -5, this will be the base of the cell
yl = 0; % bisection method
yr = 0.5;
while f_shape(yr) > -5
    yr = (1 + yr) / 2;
end
yr_init = yr;
mid = (yl+yr)/2;
fmid = f_shape(mid);
while abs(fmid - -5) > 0.001
    if fmid < -5
        yr = mid;
    else
        yl = mid;
    end
    mid = (yl+yr)/2;
    fmid = f_shape(mid);
end

total_arclength = midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048);

for i = 0:N
    target_arclength = total_arclength * i / N;
    yl = 0;
    yr = yr_init;
    mid = (yl+yr)/2;
    fmid = f_shape(mid);
    while abs(midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) - target_arclength) > total_arclength / N * 0.001
        if midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) > target_arclength
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f_shape(mid);
    end
    rv(1,N+1-i) = f_shape(mid);
    rv(2,N+1-i) = mid;
end
rv(1,N+1) = 1/a;
rv(2,N+1) = 0;

% rv(2,:) = mid:-mid/N:0;
% rv(1,:) = f(rv(2,:));
% rv(1,end) = 1/a;
rv(1,:) = rv(1,:) - min(rv(1,:));

% rv(:,2:end) = rv(:,2:end) + [-diff(rv(2,:)); diff(rv(1,:))] .* (rand(1,N) - 1/2) / 10; % add in a random perturbation in the normal direction at each point

rv_init = rv; % save the initial shape for later
% create ParabolicArc objects for the initial configuration (with dummy tension data)
initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
%[ rv, adj ] = generate_arc(N);
ii = 1;
R = 1;
[ tmp, adj ] = generate_sphere(R,N); % get the adjacency data for N patches
LUT = zeros(N,size(adj,1));
for n = 1:N
    nv = find(adj(n,:));
    nv = nv(nv>n);
    for nn = nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end
for i=N:-1:1
   index = find(LUT(i,:));
   L0Linear(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2); % linear segments
   R0Linear(i) = (rv(2,index(1))+rv(2,index(2)))/2;
   L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
   R0Parabolic(i) = initialArcs(i+1).vert(2);
end
[L0Spline R0Spline] = spline_intrinsics(rv);
L0SplineLocal = [0 cumsum(L0Spline)'];
R0SplineLocal = rv(2,:);
L0SplineLocalInit = L0SplineLocal;
R0SplineLocalInit = R0SplineLocal;
KSplineLocal = k * ones(N+1, 1);
MUSplineLocal = mu * ones(N+1, 1);

%get the bond to vertex map
nbonds = sum(adj(:))/2;
Tl = zeros(nbonds,1);% meridional tension
Tr = zeros(nbonds,1);% circumferential tension
[cl,cr] = compute_curvatures(rv);

cl_t(1).dat = cl;
cr_t(1).dat = cr;
rv_t(1).dat = rv;
LUT_t(1).dat = LUT;
adj_t(1).dat = adj;
L0_t(1).dat = L0Parabolic;
R0_t(1).dat = R0Parabolic;
Tl_t(1).dat = Tl;
Tr_t(1).dat = Tr; 

tic;

%% Simulate effects of elastic deformation
ext_verts = find(sum(adj,1)==1);

% f = zeros(num_frames, 2*(N+1));
% R0f = zeros(num_frames, N+1);
f_shape = zeros(1, 2*(N+1));
R0f = zeros(1, N+1);

%% Establish the relationship between arclength and secretion for unperturbed shape

% compute initial deformation with ODE solver + nonlinear solver
if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "spline_local"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
end

step = 1;
f(1,:) = frames(end,:); % save the established shape

for step = 2:5
    step = step + 1; % advance to next frame
    
    % Compute secretion data
    strainl = strainFramesl(end,:);
    strainr = strainFramesr(end,:);
    s = zeros(1, size(L0Spline, 2) + 1);
    gam = zeros(1, size(L0Spline, 2) + 1);
    gam_s = zeros(1, size(L0Spline, 2) + 1);
    gam_theta = zeros(1, size(L0Spline, 2) + 1);
    s(1) = 0;
    gam(1) = 1;
    gam_s(1) = strainl(end) - 1;
    gam_theta(1) = strainl(end) - 1;
    for i = 1:size(L0Spline, 1)-1
        s(i+1) = s(i) + L0Spline(N-i+1);
        gam_theta(i+1) = (rv(2,N-i+1) - rv(2,N-i+2)) / L0Spline(N-i+1) / rv(2,N-i+1) ... % dr/ds / r
            * gam_s(1:i) * L0Spline(N:-1:N-i+1); % \int_0^s \gamma_s(w) dw
        gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
        gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
    end
    velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_s(1:end-1))]); % velocity relative to the tip
    eps_s = diff(fliplr(velocity)) ./ diff(s); % combined growth rates
    eps_theta = 1 ./ fliplr(rv(2,2:end-1)) .* diff(fliplr(rv(2,2:end))) ./ diff(s) .* fliplr(velocity(1:end-1)); % TODO: double-check these indices
    
%     % update the intrinsic data
%     L0Spline = L0Spline + tstep * fliplr(gam_s)' .* L0Spline;
%     L0SplineLocal = [0 cumsum(L0Spline)'];
%     R0SplineLocal(2:N+1) = R0SplineLocal(2:N+1) + tstep * fliplr(gam_theta) .* R0SplineLocal(2:N+1);
%     L0Parabolic = L0Parabolic + tstep * fliplr(gam_s)' .* L0Parabolic;
%     R0Parabolic = R0Parabolic + tstep * fliplr(gam_s)' .* R0Parabolic;

    % Update the intrinsic data according to secretion
    [R0A R0B R0C R0D] = find_splines(L0SplineLocal, R0SplineLocal, 0, -1);
    drds = -( 3 * R0A .* L0SplineLocal(1:end-1) .^2 + 2 * R0B .* L0SplineLocal(1:end-1) + R0C );
    drds(end+1) = 1;
    
    drds_patches = -diff(rv(2,:)) ./ L0Spline';
    
    L0SplineNew = L0Spline + tstep * fliplr(fliplr(L0Spline') .* gam_s)';
%     L0SplineLocal = L0SplineLocal + tstep * L0SplineLocal .* fliplr([cumsum(fliplr(L0Spline') .* gam_s) 0]);
%     L0SplineLocal = L0SplineLocal + tstep * L0SplineLocal .* [0 cumsum(L0Spline') .* fliplr(gam_s)];
    L0SplineLocal = [0 cumsum(L0Spline)'];
%     R0SplineLocal = R0SplineLocal + tstep * drds .* R0SplineLocal .* fliplr([cumsum(fliplr(L0Spline') .* gam_s) 0]);
%     R0SplineLocal = R0SplineLocal + tstep * drds .* R0SplineLocal .* [0 cumsum(L0Spline') .* fliplr(gam_s)];
    R0SplineLocal = R0SplineLocal + tstep * drds .* fliplr([0 cumsum(fliplr(L0Spline') .* gam_s)]);
    L0Spline = L0SplineNew;
%     L0Parabolic = L0Parabolic + tstep * L0Parabolic .* fliplr(gam_s);
%     R0Parabolic = R0Parabolic + tstep * drds_patches .* fliplr(cumsum(fliplr(L0Parabolic) .* gam_s));
    
    R0f(step,:) = R0SplineLocal; % save the R0 data
    
    % Now find the new configuration using the nonlinear solver
    rv = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
    % [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
    
    % f(step,:) = frames(end,:);
    f(step,:) = reshape(rv',2*(N+1),1)';
end

s_tip = [0 cumsum(fliplr(L0Spline'))];
s_tip = s_tip(1:end-1);
[gamA gamB gamC gamD] = find_splines(s, gam, (gam(2) - gam(1)) / (s(2) - s(1)), 0);

%% Perturb the profile

a = 9;

f_shape = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
fp = @(y) pi/a .* cot(pi*y) - pi^2*y/a.*csc(pi*y).^2;

% find the y such that f(y) = -5, this will be the base of the cell
yl = 0; % bisection method
yr = 0.5;
while f_shape(yr) > -5
    yr = (1 + yr) / 2;
end
yr_init = yr;
mid = (yl+yr)/2;
fmid = f_shape(mid);
while abs(fmid - -5) > 0.001
    if fmid < -5
        yr = mid;
    else
        yl = mid;
    end
    mid = (yl+yr)/2;
    fmid = f_shape(mid);
end

total_arclength = midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048);

for i = 0:N
    target_arclength = total_arclength * i / N;
    yl = 0;
    yr = yr_init;
    mid = (yl+yr)/2;
    fmid = f_shape(mid);
    while abs(midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) - target_arclength) > total_arclength / N * 0.001
        if midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) > target_arclength
            yr = mid;
        else
            yl = mid;
        end
        mid = (yl+yr)/2;
        fmid = f_shape(mid);
    end
    rv(1,N+1-i) = f_shape(mid);
    rv(2,N+1-i) = mid;
end
rv(1,N+1) = 1/a;
rv(2,N+1) = 0;

% rv(2,:) = mid:-mid/N:0;
% rv(1,:) = f(rv(2,:));
% rv(1,end) = 1/a;
rv(1,:) = rv(1,:) - min(rv(1,:));

% create ParabolicArc objects for the initial configuration (with dummy tension data)
initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
%[ rv, adj ] = generate_arc(N);
ii = 1;
R = 1;
[ tmp, adj ] = generate_sphere(R,N); % get the adjacency data for N patches
LUT = zeros(N,size(adj,1));
for n = 1:N
    nv = find(adj(n,:));
    nv = nv(nv>n);
    for nn = nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end
for i=N:-1:1
   index = find(LUT(i,:));
   L0Linear(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2); % linear segments
   R0Linear(i) = (rv(2,index(1))+rv(2,index(2)))/2;
   L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
   R0Parabolic(i) = initialArcs(i+1).vert(2);
end
[L0Spline R0Spline] = spline_intrinsics(rv);
L0SplineLocal = [0 cumsum(L0Spline)'];
R0SplineLocal = rv(2,:);
L0SplineLocalInit = L0SplineLocal;
R0SplineLocalInit = R0SplineLocal;
KSplineLocal = k * ones(N+1, 1);
MUSplineLocal = mu * ones(N+1, 1);

% compute initial deformation with ODE solver + nonlinear solver
if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "spline_local"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
end

step = 1;
f(6,:) = frames(end,:); % save the established shape

for step = 7:num_frames
    step = step + 1; % advance to next frame
    
    strainl = strainFramesl(end,:);
    strainr = strainFramesr(end,:);
    
    % Compute secretion data
    s_current = [0 cumsum(fliplr(L0Spline'))];
    s_current = s_current(1:end-1);
    for i = 1:size(s_current, 2)-2
        s_tmp = s_current(i);
        ind = min(size(s_current, 2)-2, max(find(s_tip <= s_tmp)));
        if s_tmp > max(s_tip)
            gam(i) = 0;
        else
            gam(i) = gamA(ind) * s_tmp^3 + gamB(ind) * s_tmp^2 + gamC(ind) * s_tmp + gamD(ind);
        end
        gam_s(i) = gam(i) * strainl(end-i+1);
        gam_theta(i) = gam(i) * strainr(end-i+1);
    end
    gam(end-1:end) = [0 0];
    
    % Update the intrinsic data according to secretion
    [R0A R0B R0C R0D] = find_splines(L0SplineLocal, R0SplineLocal, 0, -1);
    drds = -( 3 * R0A .* L0SplineLocal(1:end-1) .^2 + 2 * R0B .* L0SplineLocal(1:end-1) + R0C );
    drds(end+1) = 1;
    drds_patches = -diff(rv(2,:)) ./ L0Spline';
    L0SplineNew = L0Spline + tstep * fliplr(fliplr(L0Spline') .* gam_s)';
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal = R0SplineLocal + tstep * drds .* fliplr([0 cumsum(fliplr(L0Spline') .* gam_s)]);
    L0Spline = L0SplineNew;
    R0f(step,:) = R0SplineLocal; % save the R0 data
    
    % Now find the new configuration using the nonlinear solver
    rv = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
    % [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
    
    % f(step,:) = frames(end,:);
    f(step,:) = reshape(rv',2*(N+1),1)';
end

%% Save variables, media, and cleanup

save(['growth-timelapse-hyphoid-6.8-perturbed-recovery-', num2str(a), '.mat']);

% plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);

v = VideoWriter(['media/timelapse-growth-', name, '-stiffness-', num2str(stiffness), '-perturbed-recovery-', num2str(a), '.avi']);
v.Quality = 100;
open(v);
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    % plot(p(1,:), R0f(frameNum,:), 'LineWidth', 2.0); % plot the intrinsic R0 data against the z data for this frame
    set(gca,'ColorOrderIndex',1); % reset the coloring for the rest of the plots
    slope_cutoff = -1 + 0.001; % get the spline data
    % slope_cutoff = -0.3 + 0.001; % get the spline data
    cut = 2;
    slope = (p(2,3) - p(2,1)) / (p(1,3) - p(1,1));
    while slope > slope_cutoff
        cut = cut+1;
        slope = (p(2,cut+1) - p(2,cut-1)) / (p(1,cut+1) - p(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(p(1,1:cut), p(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-p(2,cut:end), p(1,cut:end), -1/slope, 0);
    hold on; % do the plot
    for i = 1:size(horizA,2)
        spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
        plot(p(1,i):0.001:p(1,i+1), spl(p(1,i):0.001:p(1,i+1)), 'LineWidth', 2.0);
    end
    for i = 1:size(vertA,2)
        spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
        plot(spl(-p(2,cut+i-1):0.001:-p(2,cut+i)), -(-p(2,cut+i-1):0.001:-p(2,cut+i)), 'LineWidth', 2.0);
    end
    daspect([1 1 1]);
    % xlim([0 7]); ylim([0 1.5]);
    xlim([0 max(f(end,:))+0.5]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    writeVideo(v, thisFrame);
    clf(gcf);
end
close(v);

v = VideoWriter(['media/timelapse-growth-with-initial-', name, '-stiffness-', num2str(stiffness), '-perturbed-recovery-', num2str(a), '.avi']);
v.Quality = 100;
open(v);
% save the data from the first frame
pInit = reshape(f(1,:), N+1, 2)';
slope_cutoffInit = -1 + 0.001; % get the spline data
% slope_cutoffInit = -0.3 + 0.001; % get the spline data
cutInit = 2;
slopeInit = (pInit(2,3) - pInit(2,1)) / (pInit(1,3) - pInit(1,1));
while slopeInit > slope_cutoffInit
    cutInit = cutInit+1;
    slopeInit = (pInit(2,cutInit+1) - pInit(2,cutInit-1)) / (pInit(1,cutInit+1) - pInit(1,cutInit-1));
end
[horizAInit horizBInit horizCInit horizDInit] = find_splines(pInit(1,1:cutInit), pInit(2,1:cutInit), 0, slopeInit);
[vertAInit vertBInit vertCInit vertDInit] = find_splines(-pInit(2,cutInit:end), pInit(1,cutInit:end), -1/slopeInit, 0);
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    slope_cutoff = -1 + 0.001; % get the spline data
    % slope_cutoff = -0.3 + 0.001; % get the spline data
    cut = 2;
    slope = (p(2,3) - p(2,1)) / (p(1,3) - p(1,1));
    while slope > slope_cutoff
        cut = cut+1;
        slope = (p(2,cut+1) - p(2,cut-1)) / (p(1,cut+1) - p(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(p(1,1:cut), p(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-p(2,cut:end), p(1,cut:end), -1/slope, 0);
    hold on; % do the plot
    % plot for this frame's profile
    for i = 1:size(horizA,2)
        spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
        set(gca,'ColorOrderIndex',1); % reset the coloring
        plot(p(1,i):0.001:p(1,i+1), spl(p(1,i):0.001:p(1,i+1)), 'LineWidth', 2.0);
    end
    for i = 1:size(vertA,2)
        spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
        set(gca,'ColorOrderIndex',1); % reset the coloring
        plot(spl(-p(2,cut+i-1):0.001:-p(2,cut+i)), -(-p(2,cut+i-1):0.001:-p(2,cut+i)), 'LineWidth', 2.0);
    end
    % plot for the first frame's profile, displaced to pass through the
    % current tip
    for i = 1:size(horizAInit,2)
        spl = @(x) horizAInit(i) * x.^3 + horizBInit(i) * x.^2 + horizCInit(i) * x + horizDInit(i);
        set(gca,'ColorOrderIndex',2); % reset the coloring
        plot(pInit(1,i)-pInit(1,end)+p(1,end):0.001:pInit(1,i+1)-pInit(1,end)+p(1,end), spl(pInit(1,i):0.001:pInit(1,i+1)), '--', 'LineWidth', 2.0);
    end
    for i = 1:size(vertAInit,2)
        spl = @(x) vertAInit(i) * x.^3 + vertBInit(i) * x.^2 + vertCInit(i) * x + vertDInit(i);
        set(gca,'ColorOrderIndex',2); % reset the coloring
        plot(spl(-pInit(2,cutInit+i-1):0.001:-pInit(2,cutInit+i))-pInit(1,end)+p(1,end), -(-pInit(2,cutInit+i-1):0.001:-pInit(2,cutInit+i)), '--', 'LineWidth', 2.0);
    end
    daspect([1 1 1]);
    % xlim([0 7]); ylim([0 1.5]);
    xlim([0 max(f(end,:))+0.5]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    writeVideo(v, thisFrame);
    clf(gcf);
end
close(v);

function res = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun * 4,'TolX',TolX * 4,'Algorithm','levenberg-marquardt');
    % options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX);
    [X, fval, exitflag, output] = fsolve(@solver_cspline_local_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,0);
    if exitflag <= 0
        disp("fsolve failed to find solution");
    end
    res = reshape(X',N,2)'; 
end