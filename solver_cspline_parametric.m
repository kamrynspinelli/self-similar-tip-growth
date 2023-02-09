function [X_dot] = solver_cspline_local(t,X,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
global P 
%X is a row vector
%X_dot is a column vector
N = size(LUT,2); % Number of vertices
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N

r_dot = zeros(2,N); % preallocate force matrix

% interpolate the shape parametrically
[ZA ZB ZC ZD] = find_splines(L0, rv(1,:), 1, 0); % interpolate z as function of s^0
[RA RB RC RD] = find_splines(L0, rv(2,:), 0, -1); % r as function of s^0

ZAext = [ZA(1) ZA];
ZBext = [ZB(1) ZB];
ZCext = [ZC(1) ZC];
ZDext = [ZD(1) ZD];
RAext = [RA(1) RA];
RBext = [RB(1) RB];
RCext = [RC(1) RC];
RDext = [RD(1) RD];

% shape data
Zp = 3 * ZAext .* L0.^2 + 2 * ZBext .* L0 + ZCext; % dz/ds^0
Rp = 3 * RAext .* L0.^2 + 2 * RBext .* L0 + RCext; % dr/ds^0
drds = Rp ./ Zp; % dr/ds = (dr/ds^0) / (dz/ds^0)
sinalpha = ones(1,N) ./ (1 + drds.^2).^(1/2); % components of tangent vector
cosalpha = drds ./ (1 + drds.^2).^(1/2);
tvec = [sinalpha; cosalpha]; % tangent vector itself
tvec(:,end) = [0; -1];
nvec = [-cosalpha; sinalpha]; % normal vector, rotation of tangent by pi/2
nvec(:,end) = [1; 0];
Zpp = 6 * ZAext .* L0 + 2 * ZBext; % d^2r/d(s^0)^2
Rpp = 6 * RAext .* L0 + 2 * RBext; % d^2r/d(s^0)^2
ks = (Zp .* Rpp - Rp .* Zpp) ./ (Zp.^2 + Rp.^2).^(3/2); % meridional curvature
ktheta = -sinalpha ./ rv(2,:); % circumferential curvature
ktheta(end) = ks(end); % tip isotropy

% tension data
[R0A R0B R0C R0D] = find_splines(L0, R0, 0, -1);
R0Aext = [R0A(1) R0A];
R0Bext = [R0B(1) R0B];
R0Cext = [R0C(1) R0C];
R0Dext = [R0D(1) R0D];
dr0ds0 = 3 * R0Aext .* L0.^2 + 2 * R0Bext .* L0 + R0Cext; % dr^0/ds^0
Ls = (Zp.^2 + Rp.^2).^(1/2); % meridional strain
Ltheta = rv(2,:) ./ R0; % circumferential strain
Ltheta(end) = Ls(end); % tip isotropy
Ss = 1/2 * MU' .* (Ltheta.^-2 - Ls.^-2) + K' .* (Ls .* Ltheta - 1); % meridional tension
Stheta = 1/2 * MU' .* (Ls.^-2 - Ltheta.^-2) + K' .* (Ls .* Ltheta - 1); % circumferential tension
dLs = (Zp + Rp) ./ (Zp.^2 + Rp.^2); % d\lambda_s/ds
dLtheta = Ls.^-1 ./ R0 .* (Rp - dr0ds0); % d\lambda_\theta/ds
% assuming bulk and shear moduli spatially constant
dSs = MU' .* (Ls.^-3 .* dLs - Ltheta.^-2 .* dLtheta) + K' .* (Ls .* dLtheta + Ltheta .* dLs); % d\sigma_s/ds

r_dot = tvec .* dSs + nvec .* Ss .* ks + nvec .* Ss .* ktheta + P * nvec; % assemble force components

r_dot(1,1) = 0; % clamp the rear point to the r-axis
r_dot(2,end) = 0; % clamp the tip point to the z-axis

if ext_force_status == 1
    r_dot = r_dot + external_force(rv);
end

if max(any(isnan(r_dot))) || max(any(isinf(r_dot)))
    disp("force was NaN or Inf");
end

% if min(diff(rv(1,:))) < 0 || max(diff(rv(2,:))) > 0
%     disp("vertices out of order");
% end

X_dot = reshape(r_dot',2*N,1);

end

function s = signed_sqrt(a)
    s = sign(a) * sqrt(abs(a));
end