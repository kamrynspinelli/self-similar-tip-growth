function [X_dot] = solver_cspline_local(t,X,LUT,L0,R0,POIS,YOUNG,ext_verts,ext_force_status)
global P 
%X is a row vector
%X_dot is a column vector
N = size(LUT,2); % Number of vertices
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N

%% Collect all the information for the cubic spline system
% Decide which point will separate the horizontally- and
% vertically-oriented splines. This will happen when dr/dz = slope_cutoff.
% slope_cutoff = -1/2;
slope_cutoff = -1 + 0.001;
cut = 2;
slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
while slope >= slope_cutoff
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
    Ls(i) = 3 * SA(i) * S0(i)^2 + 2 * SB(i) * S0(i) + SC(i); % ds/ds^0
    d2s(i) = 6 * SA(i) * S0(i) + 2 * SB(i); % d^2(s)/d(s^0)^2
end
Ls(N) = 3 * SA(end) * S0(end)^2 + 2 * SB(end) * S0(end) + SC(end);
Ltheta(1:N-1) = rv(2,1:end-1) ./ R0(1:end-1); % need to fix last one
% Ltheta(N) = 1; % 0/0
Ltheta(N) = Ls(N); % by L'Hopital's rule
% sigmaS = YOUNG' ./ (1 - POIS') .* (Ls - 1); % WRONG
% sigmaTheta = YOUNG' ./ (1 - POIS') .* (Ltheta - 1);
sigmaS = YOUNG' ./ (2 * (1 - POIS')) .* (Ls + Ltheta - 2) + YOUNG' ./ (2 * (1 + POIS')) .* (Ls - Ltheta);
sigmaTheta = YOUNG' ./ (2 * (1 - POIS')) .* (Ls + Ltheta - 2) + YOUNG' ./ (2 * (1 + POIS')) .* (Ltheta - Ls);
% dsigmaS = YOUNG(1:end-1)' ./ (1 - POIS(1:end-1)') .* d2s .* Ls(1:end-1).^-1;
% for i = 1:N-1
%     Rp(i) = 3 * RA(i) * S(i)^2 + 2 * RB(i) * S(i) + RC(i); % dr/ds
%     R0p(i) = 3 * R0A(i) * S0(i)^2 + 2 * R0B(i) * S0(i) + R0C(i); % dr^0/ds^0
%     dsigmaS(i) = MU(i) * (Ls(i)^-4 * d2s(i) - Ltheta(i)^-3 * Rp(i) / R0(i) + Ltheta(i)^-2 * Ls(i)^-1 * R0p(i) / R0(i)) ...
%         + K(i) * (Ls(i) * Rp(i) / R0(i) - Ltheta(i) * R0p(i) / R0(i) + Ltheta(i) * Ls(i)^-1 * d2s(i));
% end
[sSA sSB sSC sSD] = find_splines(S, sigmaS, (sigmaS(2) - sigmaS(1)) / (S(2) - S(1)), (sigmaS(end) - sigmaS(end-1)) / (S(end) - S(end-1)));
for i = 1:N-1
    dsigmaS(i) = 3 * sSA(i) * S(i)^2 + 2 * sSB(i) * S(i) + sSC(i);
end
% Rp(N) = 3 * RA(end) * S(end)^2 + 2 * RB(end) * S(end) + RC(end);
% R0p(N) = 3 * R0A(end) * S0(end)^2 + 2 * R0B(end) * S0(end) + R0C(end); % dr/ds
% dsigmaS(N) = -MU(end) * (Ltheta(end)^-3 * Rp(end) / R0(end) - Ltheta(end)^-2 * Ls(end)^-1 * R0p(end) / R0(end)) ...
%     + K(end) * (Ls(end) * Rp(end) / R0(end) - Ltheta(end) * R0p(end) / R0(end));
% The last entry of dsigmaS is not assigned. This is okay, because the \hat{t}
% term for the tip point gets thrown out anyway.

% ===== DEBUG =====
% [curv_s tmp] = compute_curvatures(rv);
% curv_s = -curv_s;
% ===== END DEBUG =====

%% Compute the force
r_dot = zeros(2,N); % preallocate matrix of force vectors
% Left of the cut point (horizontally-oriented patches)
for index = 1:cut-1
    % geometric quantities
    z = rv(1,index); % coordinates of the point
    r = rv(2,index);
    r_p = 3 * horizA(index) * z^2 + 2 * horizB(index) * z + horizC(index); % r'(z)
    r_pp = 6 * horizA(index) * z + 2 * horizB(index); % r''(z)
    sin_alpha = 1 / sqrt(1 + r_p^2); % local angle 
    cos_alpha = r_p / sqrt(1 + r_p^2);
    t = [sin_alpha; cos_alpha]; % tangent and normal vectors at the marker point
    n = [-cos_alpha; sin_alpha];
    curv_s(index) = r_pp / (1 + r_p^2)^(3/2); % \kappa_s
    % mechanical quantities
    sS = sigmaS(index);
    sTheta = sigmaTheta(index);
    sS_p = dsigmaS(index);
    % add up the force components
    r_dot(:,index) = r_dot(:,index) + sS_p * t + sS * curv_s(index) * n ...
        + [0; (sS - sTheta) / r] ...
        - sS * sin_alpha / r * n ...
        + P * n;
end

% Cut point
% geometric quantities - use the left-hand patch for convenience
z = rv(1,cut); % coordinates of the point
r = rv(2,cut);
r_p = 3 * horizA(cut-1) * z^2 + 2 * horizB(cut-1) * z + horizC(cut-1); % r'(z)
r_pp = 6 * horizA(cut-1) * z + 2 * horizB(cut-1); % r''(z)
sin_alpha = 1 / sqrt(1 + r_p^2); % local angle 
cos_alpha = r_p / sqrt(1 + r_p^2);
t = [sin_alpha; cos_alpha]; % tangent and normal vectors at the marker point
n = [-cos_alpha; sin_alpha];
curv_s(cut) = r_pp / (1 + r_p^2)^(3/2); % \kappa_s
% mechanical quantities
sS = sigmaS(cut);
sTheta = sigmaTheta(cut);
sS_p = dsigmaS(cut);
% add up the force components
r_dot(:,cut) = r_dot(:,cut) + sS_p * t + sS * curv_s(cut) * n ...
    + [0; (sS - sTheta) / r] ...
    - sS * sin_alpha / r * n ...
    + P * n;

% Right of the cut point (vertically-oriented patches)
for index = cut+1:N-1
    % geometric quantities
    z = rv(1,index); % coordinates of the point
    r = rv(2,index);
    z_p = 3 * vertA(index-cut) * (-r)^2 + 2 * vertB(index-cut) * (-r) + vertC(index-cut); % r'(z)
    z_pp = 6 * vertA(index-cut) * (-r) + 2 * vertB(index-cut); % r''(z)
    sin_alpha = z_p / sqrt(1 + z_p^2); % local angle 
    cos_alpha = -1 / sqrt(1 + z_p^2);
    t = [sin_alpha; cos_alpha]; % tangent and normal vectors at the marker point
    n = [-cos_alpha; sin_alpha];
    curv_s(index) = z_pp / (1 + z_p^2)^(3/2); % \kappa_s
    % mechanical quantities
    sS = sigmaS(index);
    sTheta = sigmaTheta(index);
    sS_p = dsigmaS(index);
    % add up the force components
    r_dot(:,index) = r_dot(:,index) + sS_p * t + sS * curv_s(index) * n ...
        + [0; (sS - sTheta) / r] ...
        - sS * sin_alpha / r * n ...
        + P * n;
end

% Tip point
% geometric quantities
z = rv(1,N); % coordinates of the point
r = rv(2,N);
z_p = 3 * vertA(N-cut) * (-r)^2 + 2 * vertB(N-cut) * (-r) + vertC(N-cut); % r'(z)
z_pp = 6 * vertA(N-cut) * (-r) + 2 * vertB(N-cut); % r''(z)
sin_alpha = 0; % local angle 
cos_alpha = -1;
t = [sin_alpha; cos_alpha]; % tangent and normal vectors at the marker point
n = [-cos_alpha; sin_alpha];
curv_s(N) = z_pp / (1 + z_p^2)^(3/2); % \kappa_s
% mechanical quantities
sS = sigmaS(N);
sTheta = sigmaTheta(N);
% sS_p = dsigmaS(N);
% add up the force components
r_dot(:,N) = r_dot(:,N) ... % don't bother with the \hat{t}-term, it gets thrown out anyway
    + 2 * sS * curv_s(N) * n ... % sin(\alpha) / r = \kappa_\theta = \kappa_s by tip isotropy
    + P * n;

r_dot(1,1) = 0; % clamp the rear point to the r-axis
r_dot(2,end) = 0; % clamp the tip point to the z-axis

% translate from pointwise to integral condition
% ext_arclength(N+1) = arclength(end);
% ext_arclength(2:N) = arclength(2:N);
% ext_arclength(1) = arclength(2);
% weights = mean([ext_arclength(1:end-1); ext_arclength(2:end)]);
% r_dot = [r_dot(1,:) .* weights; r_dot(2,:) .* weights];

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