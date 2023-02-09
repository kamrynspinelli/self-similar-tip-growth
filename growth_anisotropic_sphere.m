clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% simulation parameters
model = 'parabolic'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 128; %initial number of element

% Bulk Modulus
k = 1;
% Shear Modulus
mu = 1;
% Uniform material property
K = k*ones(N,1);
MU = mu*ones(N,1);

%computational parameters
Rtol = 1e-13;
Tol =  1e-3;
TolFun = 1e-16;
TolX = 1e-16;
% Inc =  1; % fast movie speed
Inc = 0.25;
% Inc = 0.05; % slow movie speed


%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated
L0Linear = ones(N,1);
R0Linear = ones(N,1);
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);
[ rv(:,N/2+1:N+1), adj ] = generate_sphere(R,N/2); % set up the endcap
rv(1,N/2+1:N+1) = rv(1,N/2+1:N+1) + 2;
rv(1,1:N/2) = 0:4/N:2-4/N; % set up tube
rv(2,1:N/2) = 1;
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
%    L0Linear(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2); % linear segments
%    R0Linear(i) = (rv(2,index(1))+rv(2,index(2)))/2;
   L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
   R0Parabolic(i) = initialArcs(i+1).vert(2);
end
%[L0Spline R0Spline] = spline_intrinsics(rv);

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

if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
end
    
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
currentArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
currentArcs = currentArcs(2:end-1);
for i = 1:size(currentArcs, 2)
    currentArclength(i) = currentArcs(i).arclength;
end

strainl = strainFramesl1(end,:);
strainr = strainFramesr1(end,:);
s = zeros(1, size(L0Parabolic, 2) + 1);
gam_s = zeros(1, size(L0Parabolic, 2) + 1);
gam_theta = zeros(1, size(L0Parabolic, 2) + 1);
s(1) = 0;
gam_s(1) = 1;
gam_theta(1) = 1;
for i = 1:size(L0Parabolic, 1)-1 % wrong - L0, R0, strain are indexed in the opposite direction. also length mismatch??
    s(i+1) = s(i) + L0Parabolic(N-i+1);
    gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
        * (strainl(N-i+1) / strainr(N-i+1)) / rv_init(2,N-i+1) ...
        - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
        * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
%     gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
%         * (strainl(N-i+1) / rv(2,N-i)) ...
%         - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
%         * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
    gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
end
plot(s, gam_s, s, gam_theta, 'LineWidth', 2.0);
daspect([1 1 1]);
title(['Anisotropic growth profile: sphere']);
xlim([0 max(s)]);
ylim([min([gam_s gam_theta]) max([gam_s gam_theta])]);
xlabel("s");
ylabel("\gamma_s(s), \gamma_\theta(s)");
exportgraphics(gcf, ['media/anisotropic_growth_sphere_k_', num2str(k), '_mu_', num2str(mu), '.png']);

close all;
etime = cputime-time

function sigmaS = sigmaS(strainS, strainR, k, mu)
    sigmaS = mu/2 * (strainR^-2 - strainS^-2) + k * (strainR * strainS - 1);
end

function sigmaTheta = sigmaTheta(strainS, strainR, k, mu)
    sigmaTheta = mu/2 * (strainS^-2 - strainR^-2) + k * (strainR * strainS - 1);
end