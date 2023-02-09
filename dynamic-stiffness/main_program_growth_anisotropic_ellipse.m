% clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

a = 1;
b = 1;
init_stiff = 12;
max_stiff = 12;
min_stiff = 8;

tstep = 0.05;
num_frames = 500;

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 64; %initial number of element

% Bulk Modulus
k = init_stiff;
% Shear Modulus
mu = init_stiff;
% Uniform material property
div1 = N/2; % border between tube and tip region
div2 = floor(2*N/3); % border between first and second tip subregions
div3 = floor(5*N/6); % burder between second and third tip subregions
K(1:div1) = max_stiff * ones(N/2, 1); % tube region
K(div1+1:div2) = min_stiff * (max_stiff / min_stiff)^(2/3) * ones(div2-div1, 1); % three subregions near the tip
K(div2+1:div3) = min_stiff * (max_stiff / min_stiff)^(1/3) * ones(div3-div2, 1);
K(div3+1:N) = min_stiff * ones(N-div3, 1);
MU = K;

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

%% Simulate effects of elastic deformation
ext_verts = find(sum(adj,1)==1);

% hold on;
% daspect([1 1 1]);

f = zeros(num_frames, 2*(N+1));

% Initial deformation with the ODE solver and nonlinear solver
% Start with uniform material to establish the shape
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

f(1,:) = frames(end,:); % save the established shape

for step = 2:num_frames 
    %% Induced growth plot
    strainl = strainFramesl(end,:);
    strainr = strainFramesr(end,:);
    s = zeros(1, size(L0Spline, 2) + 1);
    gam = zeros(1, size(L0Spline, 2) + 1);
    gam_s = zeros(1, size(L0Spline, 2) + 1);
    gam_theta = zeros(1, size(L0Spline, 2) + 1);
    s(1) = 0;
    gam(1) = 1;
    % gam_s(1) = 1;
    % gam_theta(1) = 1;
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

    % Update intrinsic data according to directional secretion rates
    L0Spline = L0Spline + L0Spline .* fliplr(gam_s)' * tstep;
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal(2:end) = R0SplineLocal(2:end) + R0SplineLocal(2:end) .* fliplr(gam_theta) * tstep;
    L0Parabolic = L0Parabolic + L0Parabolic .* fliplr(gam_s)' * tstep;
    R0Parabolic = R0Parabolic + R0Parabolic .* fliplr(gam_theta)' * tstep;
    
    % Stiffen the material according to the distribution of secretion
    K = max_stiff - (max_stiff - min_stiff) * fliplr(gam)';
    MU = max_stiff - (max_stiff - min_stiff) * fliplr(gam)';
    KSplineLocal(1) = max_stiff;
    KSplineLocal(2:end) = K;
    MUSplineLocal(1) = max_stiff;
    MUSplineLocal(2:end) = MU;
    
    % Now find the new configuration
    rv = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
    
    f(step,:) = reshape(rv',2*(N+1),1)';
end

% plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);

v = VideoWriter(['media/timelapse-growth-a-', num2str(a), '-b-', num2str(b), '-maxstiff-', num2str(max_stiff), '-minstiff-', num2str(min_stiff), '.avi']);
open(v);
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

function res = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt');
    X = fsolve(@solver_cspline_local_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,0);
    res = reshape(X',N,2)'; 
end