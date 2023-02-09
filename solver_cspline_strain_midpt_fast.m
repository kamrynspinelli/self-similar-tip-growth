function [X_dot] = solver_cspline_strain_midpt_fast(X,LUT,L0,lR0,rR0,K,MU,ext_verts,ext_force_status)
% This is the same as solver_cspline, except the circumferential strain is
% calculated by taking the midpoint of each half-patch, rather than each
% full patch. This should reduce error.

global P 
%X is a row vector
%X_dot is a column vector
N = size(LUT,2); % Number of vertices
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N

%% Collect all the information for the cubic spline system
% Decide which point will separate the horizontally- and
% vertically-oriented splines. This will happen when dr/dz = slope_cutoff.
slope_cutoff = -1;
cut = 2;
slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
while slope > slope_cutoff
    cut = cut+1;
    slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
end
[horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
[vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
% ===== DEBUG =====
% hold on;
% plot(rv(1,:), rv(2,:), 'o');
% for i = 1:size(horizA,2)
%     spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
%     plot(rv(1,i):0.001:rv(1,i+1), spl(rv(1,i):0.001:rv(1,i+1)));
% end
% for i = 1:size(vertA,2)
%     spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
%     plot(spl(-rv(2,cut+i-1):0.001:-rv(2,cut+i)), -(-rv(2,cut+i-1):0.001:-rv(2,cut+i)));
% end
% daspect([1 1 1]);
% ===== END DEBUG =====
arclength = zeros(N-1,1);
bisector = zeros(N-1,1);
rmid = zeros(N-1,1);
for i = 1:cut-1
    arclength(i) = midpt(@(z) sqrt(1 + (3 * horizA(i) * z^2 + 2 * horizB(i) * z + horizC(i))^2), rv(1,i), rv(1,i+1), 256);
    bisector(i) = midpoint_by_arclength(@(z) 3 * horizA(i) * z.^2 + 2 * horizB(i) * z + horizC(i), rv(1,i), rv(1,i+1));
    lqrtrpt(i) = midpoint_by_arclength(@(z) 3 * horizA(i) * z.^2 + 2 * horizB(i) * z + horizC(i), rv(1,i), bisector(i));
    rqrtrpt(i) = midpoint_by_arclength(@(z) 3 * horizA(i) * z.^2 + 2 * horizB(i) * z + horizC(i), bisector(i), rv(1,i+1));
    lrmid(i) = horizA(i) * lqrtrpt(i)^3 + horizB(i) * lqrtrpt(i)^2 + horizC(i) * lqrtrpt(i) + horizD(i);
    rrmid(i) = horizA(i) * rqrtrpt(i)^3 + horizB(i) * rqrtrpt(i)^2 + horizC(i) * rqrtrpt(i) + horizD(i);
end
for i = cut:N-1
    arclength(i) = midpt(@(r) sqrt(1 + (3 * vertA(i-cut+1) * r^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1))^2), -rv(2,i), -rv(2,i+1), 256);
    bisector(i) = midpoint_by_arclength(@(r) 3 * vertA(i-cut+1) * r.^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1), -rv(2,i), -rv(2,i+1));
    lqrtrpt(i) = midpoint_by_arclength(@(r) 3 * vertA(i-cut+1) * r.^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1), bisector(i), -rv(2,i+1));
    rqrtrpt(i) = midpoint_by_arclength(@(r) 3 * vertA(i-cut+1) * r.^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1), -rv(2,i), bisector(i));
    lrmid(i) = -lqrtrpt(i);
    rrmid(i) = -rqrtrpt(i);
end
Ls = arclength ./ L0;
lLtheta = lrmid ./ lR0;
rLtheta = rrmid ./ rR0;
lsigmaS = 1/2 * MU .* (lLtheta.^-2 - Ls.^-2) + K .* (Ls .* lLtheta - 1);
rsigmaS = 1/2 * MU .* (rLtheta.^-2 - Ls.^-2) + K .* (Ls .* rLtheta - 1);
lsigmaTheta = 1/2 * MU .* (Ls.^-2 - lLtheta.^-2) + K .* (Ls .* lLtheta - 1);
rsigmaTheta = 1/2 * MU .* (Ls.^-2 - rLtheta.^-2) + K .* (Ls .* rLtheta - 1);

%% Compute the force
r_dot = zeros(2,N); % preallocate matrix of force vectors
% Rear point
rA = horizA(1);
rB = horizB(1);
rC = horizC(1);
rD = horizD(1);
rSigmaS = rsigmaS(1);
rSigmaTheta = rsigmaTheta(1);
rZmid = bisector(1);
% Tvec2D
r_dot(:,1) = r_dot(:,1) + [0; ...
    2 * rSigmaS * (3 * rA * rZmid^2 + 2 * rB * rZmid + rC) / sqrt(1 + (3 * rA * rZmid^2 + 2 * rB * rZmid + rC)^2)];
% Tvec3DrR
r_dot(:,1) = r_dot(:,1) + 2 * (rSigmaS - rSigmaTheta) * [0; ...
    midpt(@(z) sqrt(1 + (3 * rA * z^2 + 2 * rB * z + rC)^2) / (rA * z^3 + rB * z^2 + rC * z + rD), 0, rZmid, 16)];
% Tvec3DnR
r_dot(:,1) = r_dot(:,1) - 2 *  rSigmaS * [0; ...
    midpt(@(z) 1 / ((rA * z^3 + rB * z^2 + rC * z + rD) ...
    * sqrt(1 + (3 * rA * z^2 + 2 * rB * z + rC)^2)), 0, rZmid, 16)];
% Pvec
r_dot(:,1) = r_dot(:,1) + 2 * P * [0; rZmid - 0];

% Points adjacent to only horizontally-oriented patches
for index = 2:cut-1
    lA = horizA(index-1);
    lB = horizB(index-1);
    lC = horizC(index-1);
    lD = horizD(index-1);
    rA = horizA(index);
    rB = horizB(index);
    rC = horizC(index);
    rD = horizD(index);
    lSigmaS = lsigmaS(index-1);
    rSigmaS = rsigmaS(index);
    lSigmaTheta = lsigmaTheta(index-1);
    rSigmaTheta = rsigmaTheta(index);
    lZmid = bisector(index-1);
    rZmid = bisector(index);
    % Tvec2D
    r_dot(:,index) = r_dot(:,index) + [rSigmaS * 1 / sqrt(1 + (3 * rA * rZmid^2 + 2 * rB * rZmid + rC)^2) ...
        - lSigmaS * 1 / sqrt(1 + (3 * lA * lZmid^2 + 2 * lB * lZmid + lC)^2); ...
        rSigmaS * (3 * rA * rZmid^2 + 2 * rB * rZmid + rC) / sqrt(1 + (3 * rA * rZmid^2 + 2 * rB * rZmid + rC)^2) ...
        - lSigmaS * (3 * lA * lZmid^2 + 2 * lB * lZmid + lC) / sqrt(1 + (3 * lA * lZmid^2 + 2 * lB * lZmid + lC)^2)];
    % Tvec3DrR
    r_dot(:,index) = r_dot(:,index) + (lSigmaS - lSigmaTheta) * [0; ...
        midpt(@(z) sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2) / (lA * z^3 + lB * z^2 + lC * z + lD), lZmid, rv(1,index), 16)] ...
        + (rSigmaS - rSigmaTheta) * [0; ...
        midpt(@(z) sqrt(1 + (3 * rA * z^2 + 2 * rB * z + rC)^2) / (rA * z^3 + rB * z^2 + rC * z + rD), rv(1,index), rZmid, 16)];
    % Tvec3DnR
    r_dot(:,index) = r_dot(:,index) - lSigmaS * [midpt(@(z) -(3 * lA * z^2 + 2 * lB * z + lC) ...
        / ((lA * z^3 + lB * z^2 + lC * z + lD) * sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2)), lZmid, rv(1,index), 16); ...
        midpt(@(z) 1 / ((lA * z^3 + lB * z^2 + lC * z + lD) ...
        * sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2)), lZmid, rv(1,index), 16)] ...
        - rSigmaS * [midpt(@(z) -(3 * rA * z^2 + 2 * rB * z + rC) ...
        / ((rA * z^3 + rB * z^2 + rC * z + rD) * sqrt(1 + (3 * rA * z^2 + 2 * rB * z + rC)^2)), rv(1,index), rZmid, 16); ...
        midpt(@(z) 1 / ((rA * z^3 + rB * z^2 + rC * z + rD) ...
        * sqrt(1 + (3 * rA * z^2 + 2 * rB * z + rC)^2)), rv(1,index), rZmid, 16)];
    % Pvec
    r_dot(:,index) = r_dot(:,index) + P * [(lA * lZmid^3 + lB * lZmid^2 + lC * lZmid + lD) - rv(2,index); ...
        rv(1,index) - lZmid] + P * [rv(2,index) - (rA * rZmid^3 + rB * rZmid^2 + rC * rZmid + rD); ...
        rZmid - rv(1,index)];
end

% Cut point
lA = horizA(end);
lB = horizB(end);
lC = horizC(end);
lD = horizD(end);
rA = vertA(1);
rB = vertB(1);
rC = vertC(1);
rD = vertD(1);
lSigmaS = lsigmaS(cut-1);
rSigmaS = rsigmaS(cut);
lSigmaTheta = lsigmaTheta(cut-1);
rSigmaTheta = rsigmaTheta(cut);
lZmid = bisector(cut-1);
rRmid = bisector(cut);
% Tvec2D - wrong
r_dot(:,cut) = r_dot(:,cut) + [rSigmaS * (3 * rA * rRmid^2 + 2 * rB * rRmid + rC) / sqrt(1 + (3 * rA * rRmid^2 + 2 * rB * rRmid + rC)^2) ...
    - lSigmaS * 1 / sqrt(1 + (3 * lA * lZmid^2 + 2 * lB * lZmid + lC)^2); ...
    rSigmaS * -1 / sqrt(1 + (3 * rA * rRmid^2 + 2 * rB * rRmid + rC)^2) ...
    - lSigmaS * (3 * lA * lZmid^2 + 2 * lB * lZmid + lC) / sqrt(1 + (3 * lA * lZmid^2 + 2 * lB * lZmid + lC)^2)];
% Tvec3DrR - sketchy
r_dot(:,cut) = r_dot(:,cut) + (lSigmaS - lSigmaTheta) * [0; ...
    midpt(@(z) sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2) / (lA * z^3 + lB * z^2 + lC * z + lD), lZmid, rv(1,cut), 16)] ...
    + (rSigmaS - rSigmaTheta) * [0; ...
    midpt(@(r) sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2) / r, -rv(2,cut), rRmid, 16)];
% Tvec3DnR - sketchy
r_dot(:,cut) = r_dot(:,cut) - lSigmaS * [midpt(@(z) -(3 * lA * z^2 + 2 * lB * z + lC) ...
    / ((lA * z^3 + lB * z^2 + lC * z + lD) * sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2)), lZmid, rv(1,cut), 16); ...
    midpt(@(z) 1 / ((lA * z^3 + lB * z^2 + lC * z + lD) ...
    * sqrt(1 + (3 * lA * z^2 + 2 * lB * z + lC)^2)), lZmid, rv(1,cut), 16)] ...
    + rSigmaS * [midpt(@(r) (3 * rA * r^2 + 2 * rB * r + rC) ...
    / (r * sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2)), -rv(2,cut), rRmid, 16); ...
    midpt(@(r) (3 * rA * r^2 + 2 * rB * r + rC)^2 / (r ...
    * sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2)), -rv(2,cut), rRmid, 16)];
% Pvec
r_dot(:,cut) = r_dot(:,cut) + P * [(lA * lZmid^3 + lB * lZmid^2 + lC * lZmid + lD) - rv(2,cut); ...
    rv(1,cut) - lZmid] + P * [rRmid - -rv(2,cut); ...
    (rA * rRmid^3 + rB * rRmid^2 + rC * rRmid + rD) - rv(1,cut)];


% Points adjacent to only vertically-oriented patches
for index = cut+1:N-1
    lA = vertA(index-cut);
    lB = vertB(index-cut);
    lC = vertC(index-cut);
    lD = vertD(index-cut);
    rA = vertA(index-cut+1);
    rB = vertB(index-cut+1);
    rC = vertC(index-cut+1);
    rD = vertD(index-cut+1);
    lSigmaS = lsigmaS(index-1);
    rSigmaS = rsigmaS(index);
    lSigmaTheta = lsigmaTheta(index-1);
    rSigmaTheta = rsigmaTheta(index);
    lRmid = bisector(index-1);
    rRmid = bisector(index);
    % Tvec2D - WRONG
    r_dot(:,index) = r_dot(:,index) + [rSigmaS * (3 * rA * rRmid^2 + 2 * rB * rRmid + rC) / sqrt(1 + (3 * rA * rRmid^2 + 2 * rB * rRmid + rC)^2) ...
        - lSigmaS * (3 * lA * lRmid^2 + 2 * lB * lRmid + lC) / sqrt(1 + (3 * lA * lRmid^2 + 2 * lB * lRmid + lC)^2); ...
        rSigmaS * -1 / sqrt(1 + (3 * rA * rRmid^2 + 2 * rB * rRmid + rC)^2) ...
        - lSigmaS * -1 / sqrt(1 + (3 * lA * lRmid^2 + 2 * lB * lRmid + lC)^2)];
    % Tvec3DrR
    r_dot(:,index) = r_dot(:,index) + (lSigmaS - lSigmaTheta) * [0; ...
        midpt(@(r) sqrt(1 + (3 * lA * r^2 + 2 * lB * r + lC)^2) / r, lRmid, -rv(2,index), 16)] ...
        + (rSigmaS - rSigmaTheta) * [0; ...
        midpt(@(r) sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2) / r, -rv(2,index), rRmid, 16)];
    % Tvec3DnR - WRONG
    r_dot(:,index) = r_dot(:,index) + lSigmaS * [midpt(@(r) (3 * lA * r^2 + 2 * lB * r + lC) ...
        / (r * sqrt(1 + (3 * lA * r^2 + 2 * lB * r + lC)^2)), lRmid, -rv(2,index), 16); ...
        midpt(@(r) (3 * lA * r^2 + 2 * lB * r + lC)^2 / (r ...
        * sqrt(1 + (3 * lA * r^2 + 2 * lB * r + lC)^2)), lRmid, -rv(2,index), 16)] ...
        + rSigmaS * [midpt(@(r) (3 * rA * r^2 + 2 * rB * r + rC) ...
        / (r * sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2)), -rv(2,index), rRmid, 16); ...
        midpt(@(r) (3 * rA * r^2 + 2 * rB * r + rC)^2 / (r ...
        * sqrt(1 + (3 * rA * r^2 + 2 * rB * r + rC)^2)), -rv(2,index), rRmid, 16)];
    % Pvec - WRONG
    r_dot(:,index) = r_dot(:,index) + P * [-rv(2,index) - lRmid; ...
        -rv(1,index) - (rA * lRmid^3 + rB * lRmid^2 + rC * lRmid + rD)] + P * [rRmid - -rv(2,index); ...
        (rA * rRmid^3 + rB * rRmid^2 + rC * rRmid + rD) - -rv(1,index)];
end

% Tip point
lA = vertA(end);
lB = vertB(end);
lC = vertC(end);
lD = vertD(end);
lSigmaS = lsigmaS(end);
lSigmaTheta = lsigmaTheta(end);
lRmid = bisector(end);
% Tvec2D
r_dot(:,N) = r_dot(:,N) ...
    + [- 2 * lSigmaS * (3 * lA * lRmid^2 + 2 * lB * lRmid + lC) / sqrt(1 + (3 * lA * lRmid^2 + 2 * lB * lRmid + lC)^2); 0];
% Tvec3DrR
% no contribution in the z-direction
% Tvec3DnR
r_dot(:,N) = r_dot(:,N) + 2 * lSigmaS * [midpt(@(r) (3 * lA * r^2 + 2 * lB * r + lC) ...
    / (r * sqrt(1 + (3 * lA * r^2 + 2 * lB * r + lC)^2)), lRmid, -rv(2,N), 16); 0];
% Pvec
r_dot(:,N) = r_dot(:,N) + 2 * P * [-rv(2,N) - lRmid; 0];

r_dot(2,end) = 0; % clamp the last point to the z-axis
if ext_force_status == 1
    r_dot = r_dot + external_force(rv);
end

if max(any(isnan(r_dot))) || max(any(isinf(r_dot)))
    disp("force was NaN or Inf");
end

X_dot = reshape(r_dot',2*N,1);

end

function s = signed_sqrt(a)
    s = sign(a) * sqrt(abs(a));
end