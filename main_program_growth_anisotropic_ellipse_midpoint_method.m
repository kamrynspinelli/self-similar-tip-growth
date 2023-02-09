% clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

a = 0.7;
b = 1;
stiffness = 12;

tstep = 0.1;

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 64; %initial number of element

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
[ rv(:,N/2+1:N+1), adj ] = generate_sphere(R,N/2); % set up the endcap
rv(1,N/2+1:N+1) = a * rv(1,N/2+1:N+1) + 2;
rv(2,N/2+1:N+1) = b * rv(2,N/2+1:N+1);
rv(1,1:N/2) = 0:4/N:2-4/N; % set up tube
rv(2,1:N/2) = b;
rv_init = rv; % save the initial shape for later
% create ParabolicArc objects for the initial configuration (with dummy tension data)
initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
%[ rv, adj ] = generate_arc(N);
LUT = zeros(N,size(adj,1));
ii = 1;
[ tmp, adj ] = generate_sphere(R,N); % get the adjacency data for N patches
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

% hold on;
% daspect([1 1 1]);

f = zeros(120, 2*(N+1));
R0f = zeros(120, N+1);

for step = 1:120
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
    
%     if mod(step-1, 10) == 0
%         plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);
%     end
    f(step,:) = frames(end,:);

    Eg(1) = error;
    %record initial condition
    [cl,cr] = compute_curvatures(rv);
    cl_t(2).dat = cl;
    cr_t(2).dat = cr;
    rv_t(2).dat = rv;
    LUT_t(2).dat = LUT;
    adj_t(2).dat = adj;
    L0_t(2).dat = L0Parabolic;
    R0_t(2).dat = R0Parabolic;
    Tl_t(2).dat = Tl;
    Tr_t(2).dat = Tr;

    toc;

    %% Induced growth plot
    % currentArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
    % currentArcs = currentArcs(2:end-1);
    % for i = 1:size(currentArcs, 2)
    %     currentArclength(i) = currentArcs(i).arclength;
    % end

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
        gam_theta(i+1) = (rv(2,N-i+1) - rv(2,N-i+2)) / L0Spline(N-i+1) / rv(2,N-i+1) ...
            * gam_s(1:i) * L0Spline(N:-1:N-i+1);
        gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
        gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
    end
    velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_s(1:end-1))]); % velocity relative to the tip
    eps_s = diff(fliplr(velocity)) ./ diff(s); % combined growth rates
    eps_theta = 1 ./ fliplr(rv(2,2:end-1)) .* diff(fliplr(rv(2,2:end))) ./ diff(s) .* fliplr(velocity(1:end-1)); % TODO: double-check these indices
    
    % Midpoint method
    % compute f(t+h/2)
    L0SplineTmp = L0Spline + L0Spline .* fliplr(gam_s)' * tstep / 2;
    L0SplineLocalTmp = [0 cumsum(L0Spline)'];
    R0SplineLocalTmp(1) = R0SplineLocal(1);
    R0SplineLocalTmp(2:N+1) = R0SplineLocal(2:N+1) + R0SplineLocal(2:N+1) .* fliplr(gam_theta) * tstep / 2;
    L0ParabolicTmp = L0Parabolic + L0Parabolic .* fliplr(gam_s)' * tstep / 2;
    R0ParabolicTmp = R0Parabolic + R0Parabolic .* fliplr(gam_theta)' * tstep / 2;
    [ TlTmp,TrTmp,rvTmp,errorTmp,framesTmp,strainFrameslTmp,strainFramesrTmp ] = equilibrium_cspline_local(rv,LUT,L0SplineLocalTmp,R0SplineLocalTmp,KSplineLocal,MUSplineLocal,ext_verts,0,L0ParabolicTmp,R0ParabolicTmp,K,MU);
    % compute the secretion data in this intermediate state
    strainlTmp = strainFrameslTmp(end,:);
    strainrTmp = strainFramesrTmp(end,:);
    sTmp = zeros(1, size(L0SplineTmp, 2) + 1);
    gamTmp = zeros(1, size(L0SplineTmp, 2) + 1);
    gam_sTmp = zeros(1, size(L0SplineTmp, 2) + 1);
    gam_thetaTmp = zeros(1, size(L0SplineTmp, 2) + 1);
    sTmp(1) = 0;
    gamTmp(1) = 1;
    gam_sTmp(1) = strainlTmp(end) - 1;
    gam_thetaTmp(1) = strainlTmp(end) - 1;
    for i = 1:size(L0SplineTmp, 1)-1
        sTmp(i+1) = sTmp(i) + L0SplineTmp(N-i+1);
        gam_thetaTmp(i+1) = (rvTmp(2,N-i+1) - rvTmp(2,N-i+2)) / L0SplineTmp(N-i+1) / rvTmp(2,N-i+1) ...
            * gam_sTmp(1:i) * L0SplineTmp(N:-1:N-i+1);
        gam_sTmp(i+1) = gam_thetaTmp(i+1) * (strainlTmp(N-i+1) - 1) / (strainrTmp(N-i+1) - 1);
        gamTmp(i+1) = gam_thetaTmp(i+1) / (strainrTmp(N-i+1) - 1);
    end
    % update the intrinsic data according to midpoint method
    L0Spline = L0Spline + tstep * fliplr(gam_sTmp)' .* L0Spline;
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal(2:N+1) = R0SplineLocal(2:N+1) + tstep * fliplr(gam_thetaTmp) .* R0SplineLocal(2:N+1);
    L0Parabolic = L0Parabolic + tstep * fliplr(gam_sTmp)' .* L0Parabolic;
    R0Parabolic = R0Parabolic + tstep * fliplr(gam_sTmp)' .* R0Parabolic;
    
    R0f(step,:) = R0SplineLocal; % save the R0 data
end

% plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);

v = VideoWriter(['media/timelapse-growth-a-', num2str(a), '-b-', num2str(b), '-stiffness-', num2str(stiffness), '.avi']);
v.Quality = 100;
open(v);
for frameNum=1:size(f,1)
    plot(rv(1,:), R0f(frameNum,:), 'LineWidth', 2.0); % plot the intrinsic R0 data
    set(gca,'ColorOrderIndex',1); % reset the coloring for the rest of the plots
    
    p = reshape(f(frameNum,:), N+1, 2)';
    slope_cutoff = -2 + 0.001; % get the spline data
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
    xlim([0 4]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

v = VideoWriter(['media/timelapse-growth-with-initial-a-', num2str(a), '-b-', num2str(b), '-stiffness-', num2str(stiffness), '.avi']);
v.Quality = 100;
open(v);
% save the data from the first frame
pInit = reshape(f(1,:), N+1, 2)';
slope_cutoffInit = -1 + 0.001; % get the spline data
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
    xlim([0 4]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);