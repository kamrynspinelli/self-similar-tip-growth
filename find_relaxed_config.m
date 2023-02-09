function rv = find_relaxed_config(L0, R0, model)
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% make sure that L0 and R0 have the right orientation
if size(L0, 1) == 1
    L0 = L0';
end
if size(R0, 1) == 1
    R0 = R0';
end

% dimensional parameters
L = 1;%characteristic length
P = 0;%characteristic pressure
N = size(L0, 1); %initial number of element

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
Inc =  1;

% external force parameters
external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
Fext = 0;
variance_factor = 1;

%% Initialize cells
% The points will initialize as a sphere - this is okay
R = 1;
[ rv, adj ] = generate_sphere(R,N);
% create ParabolicArc objects for the initial configuration (with dummy tension data)
initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
%[ rv, adj ] = generate_arc(N);
LUT = zeros(N,size(adj,1));
ii = 1;
for n = 1:N
    nv = find(adj(n,:));
    nv = nv(nv>n);
    for nn = nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end

%get the bond to vertex map
nbonds = sum(adj(:))/2;
Tl = zeros(nbonds,1);% meridional tension
Tr = zeros(nbonds,1);% circumferential tension

tic;

% Solve without growth events
ext_verts = find(sum(adj,1)==1);

if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0,R0,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_linear(rv,LUT,L0,R0,K,MU,ext_verts,0); % linear segments
end