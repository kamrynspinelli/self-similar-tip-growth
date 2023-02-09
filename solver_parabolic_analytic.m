function [X_dot] = solver_parabolic_imperative(t,X,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
global P 
%X is a row vector
%X_dot is a column vector
N = size(LUT,2); % Number of vertices
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N

%% Collect all the information for the parabolic arc system
% Extend the point and material data
ext_rv = zeros(2, N+4); % preallocate extended points matrix
% fill the first two spots with the first two points, reflected
ext_rv(1, 1:2) = -fliplr(rv(1, 2:3));
ext_rv(2, 1:2) = fliplr(rv(2, 2:3));
ext_rv(:, 3:N+2) = rv;
ext_rv(1, N+3:N+4) = fliplr(rv(1, N-2:N-1));
ext_rv(2, N+3:N+4) = -fliplr(rv(2, N-2:N-1));
ext_K = zeros(1, N+1); % preallocate extended patch parameter array
ext_K(1) = K(1); % fill in
ext_K(2:N) = K(:);
ext_K(N+1) = K(N-1);
ext_mu = zeros(1, N+1);
ext_mu(1) = MU(1);
ext_mu(2:N) = MU(:);
ext_mu(N+1) = MU(N-1);
ext_charS = zeros(1, N+1);
ext_charS(1) = L0(1);
ext_charS(2:N) = L0(:);
ext_charS(N+1) = L0(N-1);
ext_charTheta = zeros(1, N+1);
ext_charTheta(1) = R0(1);
ext_charTheta(2:N) = R0(:);
ext_charTheta(N+1) = -R0(N-1); % need to negate this one

% Preallocate the arrays to hold all the data
tmin = zeros(N+1,1); % f(tmin) = p1
tmax = zeros(N+1,1); % f(tmax) = p2
tmid = zeros(N+1,1); % f(tmid) = vertex of the parabola
arclength = zeros(N+1,1); % the arclength of f between tmin and tmax
b = zeros(2,N+1); % the midpoint of the segment p1p2
L1 = zeros(N+1,1); % the parabolic stretch factor fitting p1, p2, and q1
L2 = zeros(N+1,1); % the parabolic stretch factor fitting p1, p2, and q2
L = zeros(N+1,1); % the harmonic mean of L1 and L2
ell = zeros(N+1,1); % half the length of the segment p1p2
beta = zeros(N+1,1); % the signed angle from the basis vector z to the outward 
%   normal of the segment p1p2
alpha = zeros(N+1,1); % the signed angle from the linear segment p1p2 to the basis 
%   vector r
% Tz % the z-coordinate of the vertex of the parabola
% Tr % the r-coordinate of the vertex of the parabola
% k = zeros(N+2,1); % bulk modulus
% mu = zeros(N+2,1); % shear modulus
vert = zeros(2,N+1); % the vertex of the arc (midpoint of the patch by arclength for tip-adjacent patch)
sigmaS = zeros(N+1,1); % the meridional tension on the patch
sigmaTheta = zeros(N+1,1); % the circumfrential tension on the patch
ltr = zeros(N+1,1); % boolean in i-th entry indicates whether the patch is left-to-right
for i = 1:N-1 % Save the data for each patch, except the two tip-adjacent ones
    q1 = ext_rv(:,i);
    p1 = ext_rv(:,i+1);
    p2 = ext_rv(:,i+2);
    q2 = ext_rv(:,i+3);
    b(:,i) = (p1 + p2) / 2; 
    beta(i) = atan((p2(2) - p1(2)) / (p2(1) - p1(1))) + pi/2;
    alpha(i) = pi - beta(i);
    ell(i) = norm(p1 - b(:,i), 2);
    tmax(i) = ell(i);
    tmin(i) = -ell(i);
    tmid(i) = 0;
    if p1(1) < p2(1)
        ltr(i) = 1;
    end
    q1_rot = [sin(beta(i)) * (q1(1) - b(1,i)) - cos(beta(i)) * (q1(2) - b(2,i)), ...
        cos(beta(i)) * (q1(1) - b(1,i)) + sin(beta(i)) * (q1(2) - b(2,i))];
    q2_rot = [sin(beta(i)) * (q2(1) - b(1,i)) - cos(beta(i)) * (q2(2) - b(2,i)), ...
        cos(beta(i)) * (q2(1) - b(1,i)) + sin(beta(i)) * (q2(2) - b(2,i))];
    L1(i) = signed_sqrt((ell(i)^2 - q1_rot(1)^2) / q1_rot(2));
    L2(i) = signed_sqrt((ell(i)^2 - q2_rot(1)^2) / q2_rot(2));
    L(i) = (2 * L1(i) * L2(i)) / (L1(i) + L2(i));
    if abs(L(i)) == Inf % covers degenerate and flat cases
        arclength(i) = 2 * ell(i);
    else
        arclength(i) = 1/2 * (2 * tmax(i) * sqrt(4 * tmax(i)^2 / L(i)^4 + 1) ...
            + L(i)^2 * asinh(2 * tmax(i) / L(i)^2));
    end
    vert(:,i) = [cos(beta(i)) * ell(i)^2 / (sign(L(i)) * L(i)^2) + b(1,i); ...
        sin(beta(i)) * ell(i)^2 / (sign(L(i)) * L(i)^2) + b(2,i)]; % the vertex of the parabola
    Ls = arclength(i) / ext_charS(i); % the arclength ratio ds/ds^0
    LsInv = ext_charS(i) / arclength(i); % the arclength ratio ds^0/ds
    Ltheta = vert(2,i) / ext_charTheta(i); % the circumferential ratio r/r_0
    LthetaInv = ext_charTheta(i) / vert(2,i); % the circumferential ratio r_0/r
    sigmaS(i) = 1/2 * ext_mu(i) * (LthetaInv^2 - LsInv^2) ...
        + ext_K(i) * (Ls * Ltheta - 1);
    sigmaTheta(i) = 1/2 * ext_mu(i) * (LsInv^2 - LthetaInv^2) ...
        + ext_K(i) * (Ls * Ltheta - 1);
end
% Upper tip patch
i = N;
q1 = ext_rv(:,N);
p1 = ext_rv(:,N+1);
p2 = ext_rv(:,N+3);
q2 = ext_rv(:,N+4);
tip = ext_rv(:,N+2);
b(:,i) = (p1 + p2) / 2;
beta(i) = atan((p2(2) - p1(2)) / (p2(1) - p1(1))) + pi/2;
alpha(i) = pi - beta(i);
ell(i) = norm(p1 - b(:,i), 2);
tmax(i) = 0;
tmin(i) = -ell(i);
% q1_rot = [sin(beta(i)) * (q1(1) - b(1,i)) - cos(beta(i)) * (q1(2) - b(2,i)), ...
%     cos(beta(i)) * (q1(1) - b(1,i)) + sin(beta(i)) * (q1(2) - b(2,i))];
% q2_rot = [sin(beta(i)) * (q2(1) - b(1,i)) - cos(beta(i)) * (q2(2) - b(2,i)), ...
%     cos(beta(i)) * (q2(1) - b(1,i)) + sin(beta(i)) * (q2(2) - b(2,i))];
L(i) = p1(2) / signed_sqrt(tip(1) - p1(1));
arclength(i) = 1/4 * (2 * ell(i) * sqrt(4 * ell(i)^2 / L(i)^4 + 1) ...
    + L(i)^2 * asinh(2 * ell(i) / L(i)^2));
tmid(i) = -1 * midpoint_by_arclength(@(t) 2*t / L(i)^2, 0, ell(i));
if p2(1) < p1(1) % if the patch is right-to-left
    vert(:,i) = [sin(beta(i)), cos(beta(i)); -cos(beta(i)), sin(beta(i))] ...
        * [-tmid(i) ; (-tmid(i).^2 + ell(i)^2) / (sign(L(i)) * L(i)^2)] ...
        + b(:,i);
else
    vert(:,i) = [sin(beta(i)), cos(beta(i)); -cos(beta(i)), sin(beta(i))] ...
        * [tmid(i) ; (-tmid(i).^2 + ell(i)^2) / (sign(L(i)) * L(i)^2)] ...
        + b(:,i);
    ltr(i) = 1;
end
Ls = arclength(i) / ext_charS(i); % the arclength ratio ds/ds^0
LsInv = ext_charS(i) / arclength(i); % the arclength ratio ds^0/ds
Ltheta = vert(2,i) / ext_charTheta(i); % the circumferential ratio r/r_0
LthetaInv = ext_charTheta(i) / vert(2,i); % the circumferential ratio r_0/r
sigmaS(i) = 1/2 * ext_mu(i) * (LthetaInv^2 - LsInv^2) ...
    + ext_K(i) * (Ls * Ltheta - 1);
sigmaTheta(i) = 1/2 * ext_mu(i) * (LsInv^2 - LthetaInv^2) ...
    + ext_K(i) * (Ls * Ltheta - 1);

% Lower tip patch
i = N+1;
q1 = ext_rv(:,N);
p1 = ext_rv(:,N+1);
p2 = ext_rv(:,N+3);
q2 = ext_rv(:,N+4);
tip = ext_rv(:,N+2);
b(:,i) = (p1 + p2) / 2;
beta(i) = atan((p2(2) - p1(2)) / (p2(1) - p1(1))) + pi/2;
alpha(i) = pi - beta(i);
ell(i) = norm(p1 - b(:,i), 2);
tmax(i) = ell(i);
tmin(i) = 0;
% q1_rot = [sin(beta(i)) * (q1(1) - b(1,i)) - cos(beta(i)) * (q1(2) - b(2,i)), ...
%     cos(beta(i)) * (q1(1) - b(1,i)) + sin(beta(i)) * (q1(2) - b(2,i))];
% q2_rot = [sin(beta(i)) * (q2(1) - b(1,i)) - cos(beta(i)) * (q2(2) - b(2,i)), ...
%     cos(beta(i)) * (q2(1) - b(1,i)) + sin(beta(i)) * (q2(2) - b(2,i))];
L(i) = p1(2) / signed_sqrt(tip(1) - p1(1));
arclength(i) = 1/4 * (2 * ell(i) * sqrt(4 * ell(i)^2 / L(i)^4 + 1) ...
    + L(i)^2 * asinh(2 * ell(i) / L(i)^2));
tmid(i) = midpoint_by_arclength(@(t) 2*t / L(i)^2, 0, ell(i));
if p2(1) < p1(1) % if the patch is right-to-left
    vert(:,i) = [sin(beta(i)), cos(beta(i)); -cos(beta(i)), sin(beta(i))] ...
        * [-tmid(i) ; (-tmid(i).^2 + ell(i)^2) / (sign(L(i)) * L(i)^2)] ...
        + b(:,i);
else
    vert(:,i) = [sin(beta(i)), cos(beta(i)); -cos(beta(i)), sin(beta(i))] ...
        * [tmid(i) ; (-tmid(i).^2 + ell(i)^2) / (sign(L(i)) * L(i)^2)] ...
        + b(:,i);
    ltr(i) = 1;
end
Ls = arclength(i) / ext_charS(i); % the arclength ratio ds/ds^0
LsInv = ext_charS(i) / arclength(i); % the arclength ratio ds^0/ds
Ltheta = vert(2,i) / ext_charTheta(i); % the circumferential ratio r/r_0
LthetaInv = ext_charTheta(i) / vert(2,i); % the circumferential ratio r_0/r
sigmaS(i) = 1/2 * ext_mu(i) * (LthetaInv^2 - LsInv^2) ...
    + ext_K(i) * (Ls * Ltheta - 1);
sigmaTheta(i) = 1/2 * ext_mu(i) * (LsInv^2 - LthetaInv^2) ...
    + ext_K(i) * (Ls * Ltheta - 1);

%% Compute the force on each point
r_dot = zeros(2,N); % preallocate matrix of force vectors
for i = 1:N-2 % compute the force on each point, one at a time
    % Tvec2D
    % tangent vector to the left-adjacent patch at its midpoint
    if ltr(i) == 1
        lTvec = [sin(beta(i)) + cos(beta(i)) * -2 * tmid(i) / (sign(L(i)) * L(i)^2); ...
            -cos(beta(i)) + sin(beta(i)) * -2 * tmid(i) / (sign(L(i)) * L(i)^2)];
    else
        lTvec = [-sin(beta(i)) + cos(beta(i)) * -2 * tmid(i) / (sign(L(i)) * L(i)^2); ...
            cos(beta(i)) + sin(beta(i)) * -2 * tmid(i) / (sign(L(i)) * L(i)^2)];
    end
    lTvec = lTvec / norm(lTvec); % make it a unit vector
    if ltr(i+1) == 1
        rTvec = [sin(beta(i+1)) + cos(beta(i+1)) * -2 * tmid(i+1) / (sign(L(i+1)) * L(i+1)^2); ...
            -cos(beta(i+1)) + sin(beta(i+1)) * -2 * tmid(i+1) / (sign(L(i+1)) * L(i+1)^2)];
    else
        rTvec = [-sin(beta(i+1)) + cos(beta(i+1)) * -2 * tmid(i+1) / (sign(L(i+1)) * L(i+1)^2); ...
            cos(beta(i+1)) + sin(beta(i+1)) * -2 * tmid(i+1) / (sign(L(i+1)) * L(i+1)^2)];
    end
    rTvec = rTvec / norm(rTvec);
    lSigmaS = sigmaS(i);
    rSigmaS = sigmaS(i+1);
    r_dot(:,i) = r_dot(:,i) + rSigmaS * rTvec - lSigmaS * lTvec;
    
    % Tvec3DrR
    lSigmaS = sigmaS(i); % needed patch parameters
    rSigmaS = sigmaS(i+1);
    lSigmaTheta = sigmaTheta(i);
    rSigmaTheta = sigmaTheta(i+1);
    lAlpha = alpha(i);
    rAlpha = alpha(i+1);
    lEll = ell(i);
    rEll = ell(i+1);
    lLambda = L(i);
    rLambda = L(i+1);
    lMdpt = b(:,i);
    rMdpt = b(:,i+1);
    lTmax = tmax(i);
    lTmid = tmid(i);
    rTmin = tmin(i+1);
    rTmid = tmid(i+1);
    % integrands for the left- and right-side integrals
    A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
    B = cos(lAlpha);
    C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
    D = 4 / lLambda^4;
    % these antiderivatives come from Mathematica - don't bother trying it
    % by hand
    lInt = (lSigmaS - lSigmaTheta) * (((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*lTmax) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*lTmax)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmax^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*lTmax)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmax^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))) ...
                - ((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*lTmid) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*lTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmid^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*lTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmid^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))));
    A = -sin(rAlpha) / (sign(rLambda) * rLambda^2);
    B = cos(rAlpha);
    C = rEll^2 * sin(rAlpha) / (sign(rLambda) * rLambda^2) + rMdpt(2);
    D = 4 / rLambda^4;
    rInt = (rSigmaS - rSigmaTheta) * (((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*rTmid) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*rTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*rTmid^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*rTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*rTmid^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))) ...
                - ((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*rTmin) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*rTmin)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*rTmin^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*rTmin)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*rTmin^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))));
    % compute the force
    r_dot(:,i) = r_dot(:,i) + [0; lInt+rInt];

    % Tvec3DnR
    lSigmaS = sigmaS(i); % needed patch parameters
    rSigmaS = sigmaS(i+1);
    lAlpha = alpha(i);
    rAlpha = alpha(i+1);
    lBeta = beta(i);
    rBeta = beta(i+1);
    lEll = ell(i);
    rEll = ell(i+1);
    lLambda = L(i);
    rLambda = L(i+1);
    lMdpt = b(:,i);
    rMdpt = b(:,i+1);
    lTmax = tmax(i);
    lTmid = tmid(i);
    rTmin = tmin(i+1);
    rTmid = tmid(i+1);
    % left patch
    if ltr(i) == 1
        A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
        B = cos(lAlpha);
        C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
        D = 4 / lLambda^4;
        F = sin(lAlpha);
        G = 2 * cos(lAlpha) / (sign(lLambda) * lLambda^2);
        H = cos(lAlpha);
        J = -2 * sin(lAlpha) / (sign(lLambda) * lLambda^2);
    else  
        A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
        B = cos(lAlpha);
        C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
        D = 4 / lLambda^4;
        F = -sin(lAlpha);
        G = 2 * cos(lAlpha) / (sign(lLambda) * lLambda^2);
        H = -cos(lAlpha);
        J = -2 * sin(lAlpha) / (sign(lLambda) * lLambda^2);
    end
    lAntiz = @(t) -lSigmaS * -((G*J*asinh(sqrt(D)*t))/(A*sqrt(D)) + (((B - sqrt(B^2 - 4*A*C))*(A*G*H + A*F*J - B*G*J) - 2*A*(A*F*H - C*G*J))* ...
       atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
      (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + ...
     ((2*A*(A*F*H - C*G*J) + (B + sqrt(B^2 - 4*A*C))*(B*G*J - A*(G*H + F*J)))* ...
       atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
      (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)));
    lAntir = @(t) -lSigmaS * (((2*G^2*asinh(sqrt(D)*t))/sqrt(D) + (sqrt(2)*(2*A^2*F^2 + B*(B - sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))*log(-B + sqrt(B^2 - 4*A*C) - 2*A*t))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) - (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(B + sqrt(B^2 - 4*A*C) + 2*A*t))/(sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)) + ...
      (sqrt(2)*(-2*A^2*F^2 + B*(-B + sqrt(B^2 - 4*A*C))*G^2 + 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(-2*A*sqrt(B^2 - 4*A*C) - B^2*D*t + 4*A*C*D*t + B*sqrt(B^2 - 4*A*C)*D*t - sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(-(B^2*D*t) - B*sqrt(B^2 - 4*A*C)*D*t + 2*A*(sqrt(B^2 - 4*A*C) + 2*C*D*t) + sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)))/(2*A));
    lInt = [lAntiz(lTmax) - lAntiz(lTmid); lAntir(lTmax) - lAntir(lTmid)];
    % right patch
    if ltr(i+1) == 1
        A = -sin(rAlpha) / (sign(rLambda) * rLambda^2);
        B = cos(rAlpha);
        C = rEll^2 * sin(rAlpha) / (sign(rLambda) * rLambda^2) + rMdpt(2);
        D = 4 / rLambda^4;
        F = sin(rAlpha);
        G = 2 * cos(rAlpha) / (sign(rLambda) * rLambda^2);
        H = cos(rAlpha);
        J = -2 * sin(rAlpha) / (sign(rLambda) * rLambda^2);
    else  
        A = -sin(rAlpha) / (sign(rLambda) * rLambda^2);
        B = cos(rAlpha);
        C = rEll^2 * sin(rAlpha) / (sign(rLambda) * rLambda^2) + rMdpt(2);
        D = 4 / rLambda^4;
        F = -sin(rAlpha);
        G = 2 * cos(rAlpha) / (sign(rLambda) * rLambda^2);
        H = -cos(rAlpha);
        J = -2 * sin(rAlpha) / (sign(rLambda) * rLambda^2);
    end
    rAntiz = @(t) -rSigmaS * -((G*J*asinh(sqrt(D)*t))/(A*sqrt(D)) + (((B - sqrt(B^2 - 4*A*C))*(A*G*H + A*F*J - B*G*J) - 2*A*(A*F*H - C*G*J))* ...
       atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
      (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + ...
     ((2*A*(A*F*H - C*G*J) + (B + sqrt(B^2 - 4*A*C))*(B*G*J - A*(G*H + F*J)))* ...
       atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
      (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)));
    rAntir = @(t) -rSigmaS * (((2*G^2*asinh(sqrt(D)*t))/sqrt(D) + (sqrt(2)*(2*A^2*F^2 + B*(B - sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))*log(-B + sqrt(B^2 - 4*A*C) - 2*A*t))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) - (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(B + sqrt(B^2 - 4*A*C) + 2*A*t))/(sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)) + ...
      (sqrt(2)*(-2*A^2*F^2 + B*(-B + sqrt(B^2 - 4*A*C))*G^2 + 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(-2*A*sqrt(B^2 - 4*A*C) - B^2*D*t + 4*A*C*D*t + B*sqrt(B^2 - 4*A*C)*D*t - sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
        log(-(B^2*D*t) - B*sqrt(B^2 - 4*A*C)*D*t + 2*A*(sqrt(B^2 - 4*A*C) + 2*C*D*t) + sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
       (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)))/(2*A));
    rInt = [rAntiz(rTmid) - rAntiz(rTmin); rAntir(rTmid) - rAntir(rTmin)];
    
    % compute the force
    r_dot(:,i) = r_dot(:,i) + lInt + rInt;

    % Pvec
    lBeta = beta(i);
    rBeta = beta(i+1);
    lLambda = L(i);
    rLambda = L(i+1);
    lTmax = tmax(i);
    lTmid = tmid(i);
    rTmin = tmin(i+1);
    rTmid = tmid(i+1);
    if ltr(i) == 1 % if the left patch goes left-to-right
        lIntegrand = @(t) [cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2); ...
            sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2)] ...
            / sqrt((cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2 ...
            + (sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2) ...
            * sqrt(1 + 4*t^2 / (lLambda^4)); % the integrand for the integral over the left half-patch
    else % right-to-left
        lIntegrand = @(t) [-cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2); ...
            -sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2)] ...
            / sqrt((-cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2 ...
            + (-sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2) ...
            * sqrt(1 + 4*t^2 / (lLambda^4)); % the integrand for the integral over the left half-patch
    end
    if ltr(i+1) == 1  % if the right patch goes left-to-right
        rIntegrand = @(t) [cos(rBeta) + sin(rBeta) * 2 * t / (sign(rLambda) * rLambda^2); ...
            sin(rBeta) - cos(rBeta) * 2 * t / (sign(rLambda) * rLambda^2)] ...
            / sqrt((cos(rBeta) + sin(rBeta) * 2 * t / (sign(rLambda) * rLambda^2))^2 ...
            + (sin(rBeta) - cos(rBeta) * 2 * t / (sign(rLambda) * rLambda^2))^2) ...
            * sqrt(1 + 4*t^2 / (rLambda^4)); % the integrand for the integral over the right half-patch
    else % right-to-left
        rIntegrand = @(t) [-cos(rBeta) + sin(rBeta) * 2 * t / (sign(rLambda) * rLambda^2); ...
            -sin(rBeta) - cos(rBeta) * 2 * t / (sign(rLambda) * rLambda^2)] ...
            / sqrt((-cos(rBeta) + sin(rBeta) * 2 * t / (sign(rLambda) * rLambda^2))^2 ...
            + (-sin(rBeta) - cos(rBeta) * 2 * t / (sign(rLambda) * rLambda^2))^2) ...
            * sqrt(1 + 4*t^2 / (rLambda^4)); % the integrand for the integral over the right half-patch
    end
    lInt = midpt(lIntegrand, lTmid, lTmax, 16); % the numerical integration over the left half-patch
    rInt = midpt(rIntegrand, rTmin, rTmid, 16); % the numerical integration over the left half-patch
    r_dot(:,i) = r_dot(:,i) + P * (lInt + rInt);
end

% Tip-adjacent point
i = N-1;
q1z = rv(1,end-1);
q1r = rv(2,end-1);
pz = rv(1,end);
A = -(pz - q1z) / (q1r^2);
B = pz;
rmid = midpoint_by_arclength(@(r) 2*A*r, 0, q1r);
% Tvec2D
if ltr(end-2) == 1
    lTvec = [sin(beta(end-2)) + cos(beta(end-2)) * -2 * tmid(end-2) / (sign(L(end-2)) * L(end-2)^2); ...
        -cos(beta(end-2)) + sin(beta(end-2)) * -2 * tmid(end-2) / (sign(L(end-2)) * L(end-2)^2)];
else
    lTvec = [-sin(beta(end-2)) + cos(beta(end-2)) * -2 * tmid(end-2) / (sign(L(end-2)) * L(end-2)^2); ...
        cos(beta(end-2)) + sin(beta(end-2)) * -2 * tmid(end-2) / (sign(L(end-2)) * L(end-2)^2)];
end
lTvec = lTvec / norm(lTvec); % make it a unit vector
if ltr(end-1) == 1
    rTvec = [sin(beta(end-1)) + cos(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2); ...
        -cos(beta(end-1)) + sin(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2)];
else
    rTvec = [-sin(beta(end-1)) + cos(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2); ...
        cos(beta(end-1)) + sin(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2)];
end
rTvec = rTvec / norm(rTvec);
lSigmaS = sigmaS(end-2);
rSigmaS = sigmaS(end-1);
r_dot(:,end-1) = r_dot(:,end-1) + rSigmaS * rTvec - lSigmaS * lTvec;
% Tvec3DrR
lSigmaS = sigmaS(end-2); % needed patch parameters
rSigmaS = sigmaS(end-1);
lSigmaTheta = sigmaTheta(end-2);
rSigmaTheta = sigmaTheta(end-1);
lAlpha = alpha(end-2);
lEll = ell(end-2);
lLambda = L(end-2);
lMdpt = b(:,end-2);
lTmax = tmax(end-2);
lTmid = tmid(end-2);
% integrands for the left- and right-side integrals
    A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
    B = cos(lAlpha);
    C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
    D = 4 / lLambda^4;
    % these antiderivatives come from Mathematica - don't bother trying it
    % by hand
    lPart = (lSigmaS - lSigmaTheta) * (((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*lTmax) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*lTmax)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmax^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*lTmax)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmax^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))) ...
                - ((2*sqrt(B^2 - 4*A*C)*sqrt(D)*asinh(sqrt(D)*lTmid) - sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)* ...
                atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*lTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmid^2))) +  ...
                sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*lTmid)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*lTmid^2))))/ ...
                (2*A*sqrt(B^2 - 4*A*C))));
A = -(pz - q1z) / (q1r^2);
B = pz;
rPart = (rSigmaS - rSigmaTheta) * ((sqrt(4 * A^2 * rmid^2) - atanh(sqrt(4 * A^2 * rmid^2))) ...
    - (sqrt(4 * A^2 * q1r^2) - atanh(sqrt(4 * A^2 * q1r^2))));
r_dot(:,end-1) = r_dot(:,end-1) + [0; lPart + rPart];
% Tvec3DnR
lSigmaS = sigmaS(end-2); % needed patch parameters
rSigmaS = sigmaS(end-1);
lAlpha = alpha(end-2);
lBeta = beta(end-2);
lEll = ell(end-2);
lLambda = L(end-2);
lMdpt = b(:,end-2);
lTmax = tmax(end-2);
lTmid = tmid(end-2);
% integrands for the left- and right-side integrals
if ltr(i) == 1
    A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
    B = cos(lAlpha);
    C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
    D = 4 / lLambda^4;
    F = sin(lAlpha);
    G = 2 * cos(lAlpha) / (sign(lLambda) * lLambda^2);
    H = cos(lAlpha);
    J = -2 * sin(lAlpha) / (sign(lLambda) * lLambda^2);
else  
    A = -sin(lAlpha) / (sign(lLambda) * lLambda^2);
    B = cos(lAlpha);
    C = lEll^2 * sin(lAlpha) / (sign(lLambda) * lLambda^2) + lMdpt(2);
    D = 4 / lLambda^4;
    F = -sin(lAlpha);
    G = 2 * cos(lAlpha) / (sign(lLambda) * lLambda^2);
    H = -cos(lAlpha);
    J = -2 * sin(lAlpha) / (sign(lLambda) * lLambda^2);
end
lAntiz = @(t) -lSigmaS * -((G*J*asinh(sqrt(D)*t))/(A*sqrt(D)) + (((B - sqrt(B^2 - 4*A*C))*(A*G*H + A*F*J - B*G*J) - 2*A*(A*F*H - C*G*J))* ...
   atanh((2*A + (-B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
  (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + ...
 ((2*A*(A*F*H - C*G*J) + (B + sqrt(B^2 - 4*A*C))*(B*G*J - A*(G*H + F*J)))* ...
   atanh((2*A - (B + sqrt(B^2 - 4*A*C))*D*t)/(sqrt(4*A^2 - 4*A*C*D + 2*B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2))))/ ...
  (sqrt(2)*A*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)));
lAntir = @(t) -lSigmaS * (((2*G^2*asinh(sqrt(D)*t))/sqrt(D) + (sqrt(2)*(2*A^2*F^2 + B*(B - sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))*log(-B + sqrt(B^2 - 4*A*C) - 2*A*t))/ ...
   (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) - (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
    log(B + sqrt(B^2 - 4*A*C) + 2*A*t))/(sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)) + ...
  (sqrt(2)*(-2*A^2*F^2 + B*(-B + sqrt(B^2 - 4*A*C))*G^2 + 2*A*G*(B*F - sqrt(B^2 - 4*A*C)*F + C*G))* ...
    log(-2*A*sqrt(B^2 - 4*A*C) - B^2*D*t + 4*A*C*D*t + B*sqrt(B^2 - 4*A*C)*D*t - sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
   (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B - sqrt(B^2 - 4*A*C))*D)) + (sqrt(2)*(2*A^2*F^2 + B*(B + sqrt(B^2 - 4*A*C))*G^2 - 2*A*G*(B*F + sqrt(B^2 - 4*A*C)*F + C*G))* ...
    log(-(B^2*D*t) - B*sqrt(B^2 - 4*A*C)*D*t + 2*A*(sqrt(B^2 - 4*A*C) + 2*C*D*t) + sqrt(2)*sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)*sqrt(1 + D*t^2)))/ ...
   (sqrt(B^2 - 4*A*C)*sqrt(2*A^2 - 2*A*C*D + B*(B + sqrt(B^2 - 4*A*C))*D)))/(2*A));
lPart = [lAntiz(lTmax) - lAntiz(lTmid); lAntir(lTmax) - lAntir(lTmid)];
% compute the force
A = -(pz - q1z) / (q1r^2);
B = pz;
rPart = rSigmaS * [(-asinh(2 * A * rmid)) - (-asinh(2 * A * q1r)); ...
    sqrt(4 * A^2 * rmid^2 + 1) - sqrt(4 * A^2 * q1r^2 + 1)];
r_dot(:,end-1) = r_dot(:,end-1) + lPart + rPart;
% Pvec
lBeta = beta(end-2);
lLambda = L(end-2);
lTmax = tmax(end-2);
lTmid = tmid(end-2);
if ltr(end-2) == 1 % if the left patch goes left-to-right
    lIntegrand = @(t) [cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2); ...
        sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2)] ...
        / sqrt((cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2 ...
        + (sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2) ...
        * sqrt(1 + 4*t^2 / (lLambda^4)); % the integrand for the integral over the left half-patch
else % right-to-left
    lIntegrand = @(t) [-cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2); ...
        -sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2)] ...
        / sqrt((-cos(lBeta) + sin(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2 ...
        + (-sin(lBeta) - cos(lBeta) * 2 * t / (sign(lLambda) * lLambda^2))^2) ...
        * sqrt(1 + 4*t^2 / (lLambda^4)); % the integrand for the integral over the left half-patch
end
lPart = P * midpt(lIntegrand, lTmid, lTmax, 16); % the numerical integration over the left half-patch
rPart = P * [(-rmid) - (-q1r); A * rmid^2 - A * q1r^2];
r_dot(:,end-1) = r_dot(:,end-1) + lPart + rPart;



% Tip point
sS = sigmaS(end-1);
sTheta = sigmaTheta(end-1);
lam = L(end-1);
% Tvec2D
if ltr(end-1) == 1
    lTvec = [sin(beta(end-1)) + cos(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2); ...
        -cos(beta(end-1)) + sin(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2)];
else
    lTvec = [-sin(beta(end-1)) + cos(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2); ...
        cos(beta(end-1)) + sin(beta(end-1)) * -2 * tmid(end-1) / (sign(L(end-1)) * L(end-1)^2)];
end
lTvec = lTvec / norm(lTvec); % make it a unit vector
if ltr(end) == 1
    rTvec = [sin(beta(end)) + cos(beta(end)) * -2 * tmid(end) / (sign(L(end)) * L(end)^2); ...
        -cos(beta(end)) + sin(beta(end)) * -2 * tmid(end) / (sign(L(end)) * L(end)^2)];
else
    rTvec = [-sin(beta(end)) + cos(beta(end)) * -2 * tmid(end) / (sign(L(end)) * L(end)^2); ...
        cos(beta(end)) + sin(beta(end)) * -2 * tmid(end) / (sign(L(end)) * L(end)^2)];
end
rTvec = rTvec / norm(rTvec);
lSigmaS = sigmaS(end-1);
rSigmaS = sigmaS(end);
r_dot(:,end) = r_dot(:,end) + rSigmaS * rTvec - lSigmaS * lTvec;
% Tvec3DrR: contributes no force at the tip
% Tvec3DnR
r_dot(:,end) = r_dot(:,end) + [2 * sS * asinh(2 * A * rmid); 0];
% Pvec
r_dot(:,end) = r_dot(:,end) + [2 * rmid; 0];

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