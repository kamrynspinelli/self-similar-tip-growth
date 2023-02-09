function [X_dot] = solver_cubic(t,X,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
global P 
%X is a row vector
%X_dot is a column vector
angle_0=asin(2*R0(end)/L0(end));
M = size(LUT,1);
N = size(LUT,2); % Number of vertices
ry = ones(M,2);
ry = [ry(:,1).*0 ry(:,2)]; % unit vector in the upward direction
rv = reshape(X',N,2)'; %Vertex positions in the usual format 2 by N
rb = LUT*rv'; %boundary vectors
LUTP = LUT>0;
LUTM = LUT<0;
Nb = size(ext_verts,2);% Number of boundary vertices
Mb = zeros(2,Nb); % boundary vertex tangential movement matrix
Mb(1,rv(1,ext_verts)>1e-10) = 1;
Mb(2,rv(2,ext_verts)>1e-10) = 1;
Fl = K.*((D./L0.*rm./R0)-1)+0.5*MU.*((R0.^2./rm.^2)-L0.^2./D.^2); % sigma_s
Fr = K.*((D./L0.*rm./R0)-1)+0.5*MU.*(L0.^2./D.^2-(R0.^2./rm.^2)); % sigma_theta
patch_coeffs = CubicApproximationHelper.all_cubic_coefficients(rv);
patch_midpoints = zeros(N-3); % preallocate vector of midpoints (z-values) of 
% the patches by arclength
patch_tvecs = zeros(2,N-3); % preallocate matrix whose columns are tangent
% vectors of the patches at their midpoints
for i=1:N-3 % for each patch excluding the last two,
    df = @(x) 3*patch_coeffs(1,i)*x^2 + 2*patch_coeffs(2,i)*x + patch_coeffs(3,i);
    % the derivative of the cubic approximant for the i-th patch
    patch_midpoints(i) = CubicApproximationHelper.midpoint_by_arclength(df, rv(1,i), rv(1,i+1));
    % the z-value giving the midpoint of the i-th patch by arclength
    patch_tvecs(:,i) = [1 df(patch_midpoints(i))]; % the tangent vector at 
    % the midpoints of the i-th patch, normalized to unit length
    patch_tvecs(:,i) = patch_tvecs(:,i)/norm(patch_tvecs(:,i), 2);
end

% The meridional tension in 2D
% Term in force balance equation arising from integration of
% d(\sigma_s \hat{t})/ds
rb2D = [Fl Fl].*(-patch_tvecs);
Tvec2D = rb2D';

% The 3D effect of tension anisotropy
% Term in force balance equation arising from integration of
% (\sigma_s-\sigma_\theta) \hat{r}/r
Tvec3DrR = zeros(2,N-3); % preallocate matrix to hold a vector for each marker
for i=1:N-3
    % the i-th marker point lies between the i-1 and i-th patches
    % --- left-hand patch data ---
    al = patch_coeffs(1,i-1);
    bl = patch_coeffs(2,i-1);
    cl = patch_coeffs(3,i-1);
    dl = patch_coeffs(4,i-1);
    % the coefficients of the cubic arc forming the i-1-th patch
    integrandl = @(z) sqrt(1 + (3*al*z.^2+2*bl*z+cl)^2) / (al*z.^3+bl*z.^2+cl*z+dl);
    % the integrand for the integral over the half-patch to the left of the
    % marker point
    doml = patch_midpoints(i-1):(rv(1,i)-patch_midpoints(i-1))/4:rv(1,i);
    valsl = integrandl(doml);
    % the values of the integrand at 4 evenly-spaced points over the
    % half-patch
    intl = trapz(doml, valsl);
    % the integral over the left-hand half-patch
    
    % --- right-hand patch data ---
    ar = patch_coeffs(1,i);
    br = patch_coeffs(2,i);
    cr = patch_coeffs(3,i);
    dr = patch_coeffs(4,i);
    % the coefficients of the cubic arc forming the i-th patch
    integrandr = @(z) sqrt(1 + (3*ar*z.^2+2*br*z+cr)^2) / (ar*z.^3+br*z.^2+cr*z+dr);
    % the integrand for the integral over the half-patch to the right of
    % the marker point
    domr = rv(1,i):(patch_midpoints(i)-rv(1,i))/4:patch_midpoints(i);
    valsr = integrandl(doml);
    % the values of the integrand at 4 evenly-spaced points over the
    % half-patch
    intr = trapz(domr, valsr);
    % the integral over the right-hand half-patch
    
    % --- force computation ---
    Tvec3DrR(:,i) = [0 (Fl-Fr) * (intl - intr)] % the force in the 
    % \hat{r}-direction
end

clear al; clear bl; clear cl; clear dl; clear ar; clear br; clear cr; clear dr;
clear integrandl; clear valsl; clear intl; clear integranr; clear valsr; clear intr;
% clear obsolete variables; it's convenient to reuse their names

% The 3D effect of meridional tension
% term in force balance equation resulting from integration of
% -\sigma_s * \sin(\alpha)/r \hat{n}
Tvec3DnR = zeros(2,N-3); % preallocate matrix to hold a vector for each marker
for i=1:N-3
    % the i-th marker point lies between the i-1 and i-th patches
    % --- left-hand patch data ---
    al = patch_coeffs(1,i-1);
    bl = patch_coeffs(2,i-1);
    cl = patch_coeffs(3,i-1);
    dl = patch_coeffs(4,i-1);
    % the coefficients of the cubic arc forming the i-1-th patch
    integrandl = @(z) [-(3*al*z.^2+2*bl*z+cl)/((al*z.^3+bl*z.^2+cl*z+dl)*sqrt(1+(3*al*z.^2+2*bl*z+cl)^2)), ...
       1/((al*z.^3+bl*z.^2+cl*z+dl)*sqrt(1+(3*al*z.^2+2*bl*z+cl)^2))];
    % the integrand for the integral over the half-patch to the left of the
    % marker point
    valsl = integrandl(doml);
    % the values of the integrand at 4 evenly-spaced points over the
    % half-patch
    intl = trapz(doml, valsl);
    % the integral over the left-hand half-patch
    
    % --- right-hand patch data ---
    ar = patch_coeffs(1,i);
    br = patch_coeffs(2,i);
    cr = patch_coeffs(3,i);
    dr = patch_coeffs(4,i);
    % the coefficients of the cubic arc forming the i-th patch
    integrandr = @(z) [-(3*ar*z.^2+2*br*z+cr)/((ar*z.^3+br*z.^2+cr*z+dr)*sqrt(1+(3*ar*z.^2+2*br*z+cr)^2)), ...
       1/((ar*z.^3+br*z.^2+cr*z+dr)*sqrt(1+(3*ar*z.^2+2*br*z+cr)^2))];
    % the integrand for the integral over the half-patch to the right of the
    % marker point
    valsr = integrandl(domr);
    % the values of the integrand at 4 evenly-spaced points over the
    % half-patch
    intr = trapz(domr, valsr);
    % the integral over the left-hand half-patch
    
    % --- force computation ---
    Tvec3DnR(:,i) = intl - intr; % the resulting force on this marker point
end

% The pressure in the normal direction
% pressure term of force balance equation resulting from integration of
% P \hat{n}
Pvec = zeros(2,N-3); % preallocate matrix to hold a vector for each marker
for i=1:N-3
    % the i-th marker point lies between the i-1 and i-th patches
    % --- left-hand patch data ---
    al = patch_coeffs(1,i-1);
    bl = patch_coeffs(2,i-1);
    cl = patch_coeffs(3,i-1);
    dl = patch_coeffs(4,i-1);
    % the coefficients of the cubic arc forming the i-1-th patch
    zl = patch_midpoints(i-1); % the z-coordinate of the left-hand midpoint
    
    % --- right-hand patch data ---
    ar = patch_coeffs(1,i);
    br = patch_coeffs(2,i);
    cr = patch_coeffs(3,i);
    dr = patch_coeffs(4,i);
    % the coefficients of the cubic arc forming the i-th patch
    zr = patch_midpoints(i); % the z-coordinate of the right-hand midpoint
    
    % --- force computation ---
    Pvec(1,i) = -P * ((ar*zr^3 + br*zr^2 + cr*zr + dr) - ...
        (al*zl^3 + bl*zl^2 + cl*zl + dl));
    Pvec(2,i) = P * (zr - zl);
end





end