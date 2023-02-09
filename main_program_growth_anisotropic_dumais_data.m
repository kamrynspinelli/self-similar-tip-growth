clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

name = 'root-hair-dumais';
profile_num = 10;
stiffness = 24;

tstep = 0.3;
% num_frames = 300;

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 96; %initial number of element

% Bulk Modulus
k = stiffness;
% Shear Modulus
mu = stiffness;
% Uniform material property
K = k*ones(N,1);
MU = mu*ones(N,1);

%computational parameters
Rtol = 1e-13;
Tol =  1e-3 * 0.7;
TolFun = 1e-16;
TolX = 1e-16;
% Inc =  1; % fast movie speed
% Inc = 0.25;
Inc = 0.05; % slow movie speed


%% Initialize cells
% R = 1;%the intrinsic radius, if a sphere is to be generated
load('../cell-profiles/dumais-2004/9_2dF_cell1.mat'); % load the cell 1 data Dumais sent us
profile = size(Px,2) - 0 - 5 * (profile_num - 1);
% load('../cell-profiles/dumais-2004/9_2cF_cell2.mat'); % load the cell 2 data Dumais sent us
% profile = size(Px,2) - 1 - 3 * (profile_num - 1);
X = Px(1:max(find(Px(:,profile))),profile)'; % load the x and y coordinates
Y = Py(1:max(find(Px(:,profile))),profile)'; % they trace the profile from bottom left, to tip, to top left
Y = -Py(1:max(find(Px(:,profile))),profile)'; % use the lower profile instead of the top
Y = Y - min(Y) + 1;
[AX BX CX DX] = find_splines(1:size(X,2), X, 1, -1); % interpolate the x and y coordinates with cubic splines
[AY BY CY DY] = find_splines(1:size(X,2), Y, 0, 0);
i_maxX = find(X == max(X)); % find the index of the rightmost point
if X(i_maxX-1) > X(i_maxX+1) % find the index of the next-rightmost point
    i_2ndmaxX = i_maxX-1;
else
    i_2ndmaxX = i_maxX+1;
end
% find the 'index' of the interpolated tip using quadratic formula
if abs(i_maxX - (-2 * BX(min(i_maxX, i_2ndmaxX)) - sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)))) < 1
    i_tip = (-2 * BX(min(i_maxX, i_2ndmaxX)) - sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)));
else
    i_tip = (-2 * BX(min(i_maxX, i_2ndmaxX)) + sqrt(4*BX(min(i_maxX, i_2ndmaxX))^2 - 12 * AX(min(i_maxX, i_2ndmaxX)) * CX(min(i_maxX, i_2ndmaxX)))) / (6 * AX(min(i_maxX, i_2ndmaxX)));
end
% compute the total arclength from the tip to the top left endpoint
total_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
    + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
for i = ceil(i_tip):size(X,2)-1 % add up the arclengths for all the patches from the tip to the rear
    total_arclength = total_arclength ...
        + midpt(@(t) sqrt((3 * AX(i) * t.^2 + 2 * BX(i) * t + CX(i))^2 ...
            + (3 * AY(i) * t.^2 + 2 * BY(i) * t + CY(i))^2), i, i+1, 64);
end
target_arclengths = fliplr(0:total_arclength/N:total_arclength); % the desired arclength from the tip to the i-th marker point

rv(1,1) = X(end); rv(2,1) = Y(end); % preload the tip and rear points
rv(1,N+1) = AX(floor(i_tip)) * i_tip^3 + BX(floor(i_tip)) * i_tip^2 + CX(floor(i_tip)) * i_tip + DX(floor(i_tip));
rv(2,N+1) = AY(floor(i_tip)) * i_tip^3 + BY(floor(i_tip)) * i_tip^2 + CY(floor(i_tip)) * i_tip + DY(floor(i_tip));
for i = 2:N % find the positions of the rest of the points
    targ = target_arclengths(i);
    arclength_sum = 0;
    interp = floor(i_tip); % the number of the tipmost spline
    tip_arclength = midpt(@(t) sqrt((3 * AX(floor(i_tip)) * t.^2 + 2 * BX(floor(i_tip)) * t + CX(floor(i_tip)))^2 ...
    + (3 * AY(floor(i_tip)) * t.^2 + 2 * BY(floor(i_tip)) * t + CY(floor(i_tip)))^2), i_tip, ceil(i_tip), 64);
    % find which patch has the desired arclength coordinate
    if targ > tip_arclength
        % if the target arclength is less than the arclength from the tip
        % to the adjacent marker point, need to advance through the other
        % patches
        arclength_sum = arclength_sum + tip_arclength;
        interp = floor(i_tip) + 1;
        while arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64) < targ
            arclength_sum = arclength_sum + midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
            + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, interp+1, 64);
            interp = interp + 1;
        end
        remaining_arclength = targ - arclength_sum;
        % use bisection method to find the point on this patch which gives the
        % desired arclength
        i_l = interp;
        i_r = interp + 1;
        mid = (i_l + i_r) / 2;
        guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
            + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
        while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
            if guess_arclength - remaining_arclength < 0
                i_l = mid;
            else
                i_r = mid;
            end
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), interp, mid, 64);
        end
        rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
        rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
    else % if this arclength coordinate lies on the tip interpolant, just do bisection method
        remaining_arclength = targ;
        % use bisection method to find the point on this patch which gives the
        % desired arclength
        i_l = i_tip;
        i_r = ceil(i_tip);
        mid = (i_l + i_r) / 2;
        guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
            + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
        while abs(guess_arclength - remaining_arclength) > total_arclength / N / 100
            if guess_arclength - remaining_arclength < 0
                i_l = mid;
            else
                i_r = mid;
            end
            mid = (i_l + i_r) / 2;
            guess_arclength = midpt(@(t) sqrt((3 * AX(interp) * t.^2 + 2 * BX(interp) * t + CX(interp))^2 ...
                + (3 * AY(interp) * t.^2 + 2 * BY(interp) * t + CY(interp))^2), i_tip, mid, 64);
        end
        rv(1,i) = AX(interp) * mid^3 + BX(interp) * mid^2 + CX(interp) * mid + DX(interp);
        rv(2,i) = AY(interp) * mid^3 + BY(interp) * mid^2 + CY(interp) * mid + DY(interp);
    end
end
rv(2,:) = -rv(2,:); % use the lower profile instead of the top
rv(2,:) = rv(2,:) - min(rv(2,:)); % center at (0,0)
rv(1,:) = rv(1,:) - min(rv(1,:));
rv = rv / max(rv(2,:)); % normalize to have height 1

N = size(rv,2) - 1;
% Uniform material property
K = k*ones(N,1);
MU = mu*ones(N,1);
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
f = zeros(1, 2*(N+1));
R0f = zeros(1, N+1);

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

target_length = 1.4 * rv(1,end); % grow until the tube has elongated to twice the pre-growth length

step = 1;
f(1,:) = frames(end,:); % save the established shape

% for step = 2:num_frames
while rv(1,end) < target_length
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

% plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);

v = VideoWriter(['media/timelapse-growth-', name, '-stiffness-', num2str(stiffness), '.avi']);
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
    xlim([0 target_length+1]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

v = VideoWriter(['media/timelapse-growth-with-initial-', name, '-stiffness-', num2str(stiffness), '.avi']);
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
    xlim([0 target_length+1]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

function res = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun * 4,'TolX',TolX * 4,'Algorithm','levenberg-marquardt');
    % options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX);
    X = fsolve(@solver_cspline_local_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,0);
    res = reshape(X',N,2)'; 
end