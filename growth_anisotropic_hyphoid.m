function [tip_length gam_peak ZN RN Rtip W_int ext_peak W_int_ext gam L0Spline strainl ksspline] = growth_anisotropic_hyphoid(a,stiffness)
clearvars -except a stiffness
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline', 'spline_local'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 160; %initial number of element

% Bulk Modulus
k = stiffness;
% Shear Modulus
mu = stiffness;
% Uniform material property
K = k*ones(N,1);
MU = mu*ones(N,1);

%computational parameters
Rtol = 1e-13;
Tol = 1e-3;
% Tol =  1e-3 * 0.02;
TolFun = 1e-16;
TolX = 1e-16;
% Inc =  1; % fast movie speed
% Inc = 0.25;
% Inc = 0.05; % slow movie speed
Inc = 0.01;


%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated
L0Linear = ones(N,1);
R0Linear = ones(N,1);
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);

[ tmp, adj ] = generate_sphere(R,N); % get the adjacency data for N patches

f = @(y) pi*y/a .* cot(pi*y); % the profile curve in the form x = f(y)
fp = @(y) pi/a .* cot(pi*y) - pi^2*y/a.*csc(pi*y).^2;

% find the y such that f(y) = -5, this will be the base of the cell
yl = 0; % bisection method
yr = 0.5;
while f(yr) > -5
    yr = (1 + yr) / 2;
end
yr_init = yr;
mid = (yl+yr)/2;
fmid = f(mid);
while abs(fmid - -5) > 0.001
    if fmid < -5
        yr = mid;
    else
        yl = mid;
    end
    mid = (yl+yr)/2;
    fmid = f(mid);
end

total_arclength = midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048);

% uniform point distribution
% for i = 0:N
%     target_arclength = total_arclength * i / N;
%     yl = 0;
%     yr = yr_init;
%     mid = (yl+yr)/2;
%     fmid = f(mid);
%     while abs(midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) - target_arclength) > total_arclength / N * 0.001
%         if midpt(@(y) sqrt(1 + fp(y).^2), 0, mid, 2048) > target_arclength
%             yr = mid;
%         else
%             yl = mid;
%         end
%         mid = (yl+yr)/2;
%         fmid = f(mid);
%     end
%     rv(1,N+1-i) = f(mid);
%     rv(2,N+1-i) = mid;
% end
% rv(1,N+1) = 1/a;
% rv(2,N+1) = 0;

% equal y spaced distribution
rv(2,:) = mid:-mid/N:0;
rv(1,:) = f(rv(2,:));
rv(1,end) = 1/a;

rv(1,:) = rv(1,:) - min(rv(1,:));

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

if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "spline_local"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU); % parabolic arcs
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
% currentArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
% currentArcs = currentArcs(2:end-1);
% for i = 1:size(currentArcs, 2)
%     currentArclength(i) = currentArcs(i).arclength;
% end

strainl = strainFramesl1(end,:);
strainr = strainFramesr1(end,:);
s = zeros(1, size(L0Spline, 2) + 1);
gam = zeros(1, size(L0Spline, 2) + 1);
gam_s = zeros(1, size(L0Spline, 2) + 1);
gam_theta = zeros(1, size(L0Spline, 2) + 1);
s(1) = 0; % first entry is the tip
gam(1) = 1;
% gam_s(1) = 1;
% gam_theta(1) = 1;
gam_s(1) = strainl(end) - 1;
gam_theta(1) = strainl(end) - 1;
for i = 1:size(L0Spline, 1)-1
    s(i+1) = s(i) + L0Spline(N-i+1); % this is actually s^0
%     gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
%         * (strainl(N-i+1) / strainr(N-i+1)) / rv_init(2,N-i+1) ...
%         - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
%         * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
%     gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
%         * (strainl(N-i+1) / rv(2,N-i+1)) ...
%         - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
%         * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
%     gam_theta(i+1) = (rv(2,N-i+1) - rv(2,N-i+2)) / L0Spline(N-i+1) / rv(2,N-i+1) ...
%         * gam_s(1:i) * L0Spline(N:-1:N-i+1);)
    gam_theta(i+1) = (rv_init(2,N-i+1) - rv_init(2,N-i+2)) / L0Spline(N-i+1) / rv_init(2,N-i+1) ...
        * gam_s(1:i) * L0Spline(N:-1:N-i+1);
    gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
    gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
end
s_turgid = [0 cumsum(fliplr(L0Spline(2:end)' .* strainl(3:end)))]; % post-deformation arclength
velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_s(1:end-1))]); % velocity relative to the tip; last entry is the tip
% velocity = fliplr([0 cumsum(fliplr(strainl(2:end) .* L0Spline(2:end)') .* gam_s(1:end-1))]); % velocity relative to the tip; last entry is the tip % DEBUG
eps_s = diff(fliplr(velocity)) ./ diff(s_turgid); % combined growth rates; first entry is the tip
eps_theta = 1 ./ fliplr(rv(2,2:end-1)) .* diff(fliplr(rv(2,2:end))) ./ diff(s_turgid) .* fliplr(velocity(1:end-1));
eps_se = diff(fliplr(strainl(2:end))) ./ diff(s_turgid) .* fliplr(velocity(2:end)) ./ fliplr(strainl(3:end)); % first entry is the tip
eps_thetae = diff(fliplr(strainr(2:end))) ./ diff(s_turgid) .* fliplr(velocity(2:end)) ./ fliplr(strainr(3:end)); % first entry is the tip
[ks ktheta] = compute_curvatures(rv); % estimation by finite difference for each patch
[ksspline kthetaspline] = compute_curvatures_spline(rv); % estimation by cubic splines

% infer the extensibility \Phi from the Dumais paper
nu = 1/2;
beta = 2 * nu^2 - 2 * nu + 2;
for i = 1:size(s,2) % compute the effective stress
    sigma_e(i) = sqrt(Tl(end-i+1)^2 + Tr(end-i+1)^2 - Tl(end-i+1)*Tr(end-i+1)); % first entry is the tip!
end
sigma_y = 1/2 * min(sigma_e); % yield stress
extensibility(1) = 1; % normalize to 1 at the tip
% extensibility(1) = eps_s(1) / ((Tl(end) - sigma_y) * (1 - nu) / sqrt(3 * beta - 6 * nu)); % normalize so the two computations of \dot{\epsilon_s} agree at the tip
for i = 1:size(s, 2)-1
    sigma_s(i) = Tl(end-i+1);
    sigma_theta(i) = Tr(end-i+1);
    K_dumais(i) = sqrt(beta * sigma_s(i)^2 + beta * sigma_theta(i)^2 + (beta - 6*nu) * sigma_s(i) * sigma_theta(i));
    eps_s_dumais(i) = extensibility(i) * (sigma_e(i) - sigma_y) * (sigma_s(i) - nu * sigma_theta(i)) / K_dumais(i);
    eps_theta_dumais(i) = extensibility(i) * (sigma_e(i) - sigma_y) * (sigma_theta(i) - nu * sigma_s(i)) / K_dumais(i);
    if sigma_e(i) < sigma_y % implement the indicator function
        eps_s_dumais(i) = 0;
        eps_theta_dumais(i) = 0;
    end
    I = sum(eps_s_dumais(1:i) .* fliplr(strainl(end-i+1:end) .* L0Spline(end-i+1:end)'));
    drds(i) = (rv(2,end-i) - rv(2,end-i+1)) / (L0Spline(end-i+1) * strainl(end-i+1));
    extensibility(i+1) = drds(i) / (rv(2,end-i) * (sigma_e(i) - sigma_y) * (sigma_theta(i) - nu * sigma_s(i)) / K_dumais(i)) * I;
%     I = sum(eps_s(1:i) .* fliplr(strainl(end-i+1:end) .* L0Spline(end-i+1:end)'));
%     drds = (rv(2,end-i) - rv(2,end-i+1)) / (L0Spline(end-i+1) * strainl(end-i+1));
%     extensibility(i+1) = drds / (rv(2,end-i) * (sigma_e(i) - sigma_y) * (sigma_theta - nu * sigma_s) / _dumais(i)) * I;
end

s = [0 cumsum(fliplr(L0Spline(2:end)' .* strainl(3:end)))]; % post-deformation arclength

tip_length = strainl(N/2+2:end) * L0Spline(N/2+1:end); % post-deformation arclength
gam_peak = strainl(end-find(gam == max(gam))+2:end) * L0Spline(end-find(gam == max(gam))+2:end); % arclength of peak in gamma
ext_peak = strainl(end-find(extensibility == max(extensibility))+2:end) * L0Spline(end-find(extensibility == max(extensibility))+2:end); % arclength of peak in extensibility
ZN = rv(1,end) - 2; % subtract off intrinsic tube length
RN = max(rv(2,:));
Rtip = 1 / ksspline(end);

ks_peak = strainl(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end) * L0Spline(end-find(fliplr(ksspline) == max(fliplr(ksspline)))+2:end);

% find the secretion window using integral condition
levels = fliplr(unique(gam));
s1_ind = find(gam == max(gam));
s2_ind = find(gam == max(gam));
% total_gam = sum(gam(1:N/2) .* fliplr(L0Spline(end-N/2+1:end))' .* fliplr(strainl(end-N/2+1:end)));
total_gam = sum(gam(1:N) .* fliplr(L0Spline(end-N+1:end)') .* fliplr(strainl(end-N+1:end)));
window_int = sum(gam(s1_ind:s2_ind-1) .* fliplr(L0Spline(end-s2_ind+2:end-s1_ind+1)') .* fliplr(strainl(end-s2_ind+2:end-s1_ind+1)));
l = 1;
while window_int < 0.5 * total_gam % assume that \gamma is either monotonic or has a single peak
    l = l + 1;
    p = find(gam == levels(l));
    if p < s1_ind
        s1_ind = p;
    else
        s2_ind = p;
    end
    window_int = sum(gam(s1_ind:s2_ind-1) .* fliplr(L0Spline(end-s2_ind+2:end-s1_ind+1)') .* fliplr(strainl(end-s2_ind+2:end-s1_ind+1)));
end
s1 = sum(fliplr(L0Spline(end-s1_ind+2:end)') .* fliplr(strainl(end-s1_ind+2:end)));
s2 = sum(fliplr(L0Spline(end-s2_ind+2:end)') .* fliplr(strainl(end-s2_ind+2:end)));
W_int = [s1 s2];

% find the extensibility window
levels = fliplr(unique(extensibility));
s1_ind_ext = find(extensibility == max(extensibility));
s2_ind_ext = find(extensibility == max(extensibility));
% total_gam = sum(gam(1:N/2) .* fliplr(L0Spline(end-N/2+1:end))' .* fliplr(strainl(end-N/2+1:end)));
total_ext = sum(extensibility(1:N) .* fliplr(L0Spline(end-N+1:end)') .* fliplr(strainl(end-N+1:end)));
window_int = sum(extensibility(s1_ind_ext:s2_ind_ext-1) .* fliplr(L0Spline(end-s2_ind_ext+2:end-s1_ind_ext+1)') .* fliplr(strainl(end-s2_ind_ext+2:end-s1_ind_ext+1)));
l = 1;
while window_int < 0.5 * total_ext % assume that \gamma is either monotonic or has a single peak
    l = l + 1;
    p = find(extensibility == levels(l));
    if p < s1_ind_ext
        s1_ind_ext = p;
    else
        s2_ind_ext = p;
    end
    window_int = sum(extensibility(s1_ind_ext:s2_ind_ext-1) .* fliplr(L0Spline(end-s2_ind_ext+2:end-s1_ind_ext+1)') .* fliplr(strainl(end-s2_ind_ext+2:end-s1_ind_ext+1)));
end
s1_ext = sum(fliplr(L0Spline(end-s1_ind_ext+2:end)') .* fliplr(strainl(end-s1_ind_ext+2:end)));
s2_ext = sum(fliplr(L0Spline(end-s2_ind_ext+2:end)') .* fliplr(strainl(end-s2_ind_ext+2:end)));
W_int_ext = [s1_ext s2_ext];

save(['hyphoid_a_', num2str(a), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat'], ...
    'gam', 'gam_s', 'gam_theta', 's', 'velocity', 'eps_s', 'eps_theta', 'eps_se', 'eps_thetae', 'extensibility', 'strainl', 'strainr', 'ks', 'ktheta', 'ksspline', 'kthetaspline');

% plot \gamma_s and \gamma_\theta
plot(s, gam_s, s, gam_theta, 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Anisotropic growth profile: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min([gam_s gam_theta]) max([gam_s gam_theta])]);
xlabel("s");
ylabel("\gamma_s(s), \gamma_\theta(s)");
exportgraphics(gcf, ['media/anisotropic_growth_gamma_s_theta_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot \gamma
plot(s, gam, 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Anisotropic growth profile: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min(gam) max(gam)]);
xlabel("s");
ylabel("\gamma(s)");
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, ['media/anisotropic_growth_gamma_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot \gamma and secretion window
hold on;
plot(s, gam, 'LineWidth', 2.0);
fill([s(s1_ind:s2_ind) s(s2_ind) s(s1_ind)], [gam(s1_ind:s2_ind) 0 0], [0.5 0.6 1], 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 2.0);
pbaspect([1 1 1]);
xlim([0 max(s)]);
ylim([min(gam) max(gam)]);
xlabel("s");
ylabel("\gamma(s)");
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, ['media/anisotropic_growth_secretion_window_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot \epsilon_s and \epsilon_\theta
plot(s(1:end-1), eps_s, s(1:end-1), eps_theta, 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Anisotropic growth profile: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim tight;
ylim([min([eps_s eps_theta]) max([eps_s eps_theta])]);
xlabel("s");
ylabel('$\dot{\epsilon}_s, \dot{\epsilon}_\theta$', 'Interpreter', 'latex', 'FontSize', 14);
ax = gca;
set(ax, 'fontsize', 12);
exportgraphics(gcf, ['media/anisotropic_growth_epsilon_hyphoid_a_', num2str(a) '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot strains
plot(s, fliplr(strainl(2:end)), s, fliplr(strainr(2:end)), 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Strain distributions: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min([strainl strainr]) max([strainl strainr])]);
xlabel("s");
ylabel("\lambda_s(s), \lambda_\theta(s)");

exportgraphics(gcf, ['media/anisotropic_strain_distribution_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot monotonicity condition
r = fliplr(rv_init(2,:));
r_p = diff(r) ./ fliplr(L0Spline');
r_pp = diff(r_p) ./ fliplr(L0Spline(1:end-1)');
strain_tilde = fliplr((strainl(2:end) - 1) ./ (strainr(2:end) - 1));
strain_theta_star = fliplr(strainr(2:end) - 1);
strain_theta_star_p = diff(strain_theta_star) ./ fliplr(L0Spline(1:end-1)');
Y = r_pp + (strain_tilde(1:end-1) - 1) .* r_p(1:end-1).^2 ./ r(1:end-2) - strain_theta_star_p ./ strain_theta_star(1:end-1) .* r_p(1:end-1);
hold on;
plot(s(2:end), Y, 'LineWidth', 2.0);
title("Monotonicity of \gamma")
xlabel("s")
ylabel("(r^0)'' + (\lambda - 1) ((r^0)')^2 / r^0 - \lambda_\theta' / \lambda_\theta^* (r^0)'");
p = plot(s, zeros(size(s)), '--', 'LineWidth', 2.0);
p.Color = 'black';
exportgraphics(gcf, ['media/anisotropic_monotonicity_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot all quantities wrt z
hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none', 'Padding', 'tight'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
% xlim([2 rv(1,end)+0.3]);
xlim([2 5.4]);
% pbaspect([max(rv(1,:))+0.3 1.2 1])
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma
% p4 = plot(rv(1,2:end), velocity, '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
% p4.Color = c(5,:); % plot velocity
% p4 = plot(rv(1,3:end), fliplr(eps_s) / max([eps_s eps_theta]), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \epsilon_s');
% p4.Color = c(5,:); % plot meridional strain rate
% p5 = plot(rv(1,3:end), fliplr(eps_theta) / max([eps_s eps_theta]), '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain rate \epsilon_\theta');
% p5.Color = c(6,:); % plot circumferential strain rate
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
% ylim([0 max(gam)*1.1]);
ylim([0 1.5]);
ylabel('secretion rate $\gamma$', 'Interpreter', 'latex');
% pbaspect([2 1 1]);
ax.XAxis.Visible = 'off'; % turn off this x-axis
set(get(ax,'YLabel'),'Position',get(get(ax,'YLabel'),'Position') + [0.0 0 0]);
% title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(gca, 'YTick', []);
% xlim([2 rv(1,end)+0.3]);
xlim([2 5.4]);
% pbaspect([max(rv(1,:))+0.3 1.2 1])
% line([3.2 3.3], [0 0])
yyaxis right;
ymin = min([strainl strainr]); % adjust the y bounds for the strain plot
ymax = max([strainl strainr]);
p5 = plot(rv(1,:), strainl, '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain \lambda_s');
p5.Color = c(3,:); % plot meridional strain
p6 = plot(rv(1,:), strainr, '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain \lambda_\theta');
p6.Color = c(4,:); % plot circumferential strain
% p7 = plot(rv(1,3:end), fliplr(eps_s) / max([eps_s eps_theta]) * (ymax-ymin) + ymin, '-', 'LineWidth', 2.0, 'DisplayName', '\epsilon_s');
% p7.Color = c(3,:);
% p8 = plot(rv(1,3:end), fliplr(eps_theta) / max([eps_s eps_theta]) * (ymax-ymin) + ymin, '-', 'LineWidth', 2.0, 'DisplayName', '\epsilon_\theta');
% p8.Color = c(4,:);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
% ylim([1.025 1.097]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(ax, 'YDir', 'reverse'); % flip the y-axis upside down
ylabel('stretch ratios $\lambda_s$, $\lambda_\theta$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
ax.XRuler.TickLabelGapOffset = -20;
% pbaspect([2 1 1]);
t.Units = 'inches';
width = 8 * 0.8; height = width * 1.2 / ((rv(1,end) - 2) + 0.3); % width and height for one of the panels
t.InnerPosition = [1 1 width 2*height]; % hack to set the aspect ratio because tiled layout can't do that I guess 
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot secretion and strain rates wrt z
hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
xlim([0 3.4]);
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma
% p4 = plot(rv(1,2:end), velocity, '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
% p4.Color = c(5,:); % plot velocity
set(ax, 'fontsize', 12);
ax.YColor = 'black';
ylim([0 max(gam)*1.1]);
ylabel('secretion rate $\gamma$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxis.Visible = 'off'; % turn off this x-axis
% title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(gca, 'YTick', []);
xlim([0 3.4]);
line([3.2 3.3], [0 0])
yyaxis right;
ymin = min([eps_s eps_theta]); % adjust the y bounds for the strain plot
ymax = max([eps_s eps_theta]);
p4 = plot(rv(1,3:end), fliplr(eps_s), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \epsilon_s');
p4.Color = c(3,:); % plot meridional strain rate
p5 = plot(rv(1,3:end), fliplr(eps_theta), '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain rate \epsilon_\theta');
p5.Color = c(4,:); % plot circumferential strain rate
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(ax, 'YDir', 'reverse'); % flip the y-axis upside down
ylabel('strain rates $\dot{\epsilon}_s$, $\dot{\epsilon}_\theta$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
ax.XRuler.TickLabelGapOffset = -20;
pause(0.5);
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_strain_rates_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot secretion/extensibility and strain rate \dot{\epsilon}_s decomposition wrt z
hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
xlim([0 5.6]);
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma
p4 = plot(rv(1,2:end), fliplr(extensibility), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional extensibility rate \Phi');
p4.Color = c(5,:); % plot velocity
set(ax, 'fontsize', 12);
ax.YColor = 'black';
ylim([0 max([gam extensibility])*1.1]);
ylabel('secretion rate $\gamma$, extensibility $\Phi$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxis.Visible = 'off'; % turn off this x-axis
% title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(gca, 'YTick', []);
xlim([0 5.6]);
line([3.2 3.3], [0 0])
yyaxis right;
ymin = min([eps_s eps_se gam_s]); % adjust the y bounds for the strain rate decomposition plot
ymax = max([eps_s eps_se gam_s]);
p4 = plot(rv(1,3:end), fliplr(eps_s), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \dot{\epsilon}_s');
p4.Color = c(3,:); % plot meridional strain rate
p5 = plot(rv(1,3:end), fliplr(eps_se), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of elastic deformation to meridional strain rate \dot{\epsilon}_{se}');
p5.Color = c(4,:); % plot \dot{\epsilon}_{se}
p6 = plot(rv(1,2:end), fliplr(gam_s), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of exocytosis to meridional strain rate \gamma_s');
p6.Color = c(6,:);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(ax, 'YDir', 'reverse'); % flip the y-axis upside down
ylabel('strain rate decomposition $\dot{\epsilon}_s$, $\dot{\epsilon}_{se}$, $\gamma_s$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
ax.XRuler.TickLabelGapOffset = -20;
pause(0.5);
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_strain_rate_decomposition_s_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot secretion/extensibility and strain rate \dot{\epsilon}_\theta decomposition wrt z
hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
xlim([0 5.6]);
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma
p4 = plot(rv(1,2:end), fliplr(extensibility), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional extensibility rate \Phi');
p4.Color = c(5,:); % plot velocity
set(ax, 'fontsize', 12);
ax.YColor = 'black';
ylim([0 max([gam extensibility])*1.1]);
ylabel('secretion rate $\gamma$, extensibility $\Phi$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxis.Visible = 'off'; % turn off this x-axis
% title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(gca, 'YTick', []);
xlim([0 5.6]);
line([3.2 3.3], [0 0])
yyaxis right;
ymin = min([eps_theta eps_thetae gam_theta]); % adjust the y bounds for the strain rate decomposition plot
ymax = max([eps_theta eps_thetae gam_theta]);
p4 = plot(rv(1,3:end), fliplr(eps_theta), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \dot{\epsilon}_s');
p4.Color = c(3,:); % plot meridional strain rate
p5 = plot(rv(1,3:end), fliplr(eps_thetae), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of elastic deformation to meridional strain rate \dot{\epsilon}_{se}');
p5.Color = c(4,:); % plot \dot{\epsilon}_{se}
p6 = plot(rv(1,2:end), fliplr(gam_theta), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of exocytosis to meridional strain rate \gamma_s');
p6.Color = c(6,:);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(ax, 'YDir', 'reverse'); % flip the y-axis upside down
ylabel('strain rate decomposition $\dot{\epsilon}_\theta$, $\dot{\epsilon}_{\theta e}$, $\gamma_\theta$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
ax.XRuler.TickLabelGapOffset = -20;
pause(0.5);
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_strain_rate_decomposition_theta_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot both strain rate decompositions wrt z
hold on;
t = tiledlayout(2, 1, 'TileSpacing', 'none'); % set up the layout
c = lines(8); % some colors to use
nexttile; % upper panel
hold on;
% xlabel("profile");
yyaxis left;
p1 = plot([rv_init(1,:)], [rv_init(2,:)], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot upper unturgored profile
p2 = plot([rv(1,:)], [rv(2,:)], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot upper turgored profile
ylim([0 1.2]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
xlim([2 5.6]);
yyaxis right;
ymin = min([eps_s eps_se gam_s]); % adjust the y bounds for the strain rate decomposition plot
ymax = max([eps_s eps_se gam_s]);
p4 = plot(rv(1,3:end), fliplr(eps_s), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \dot{\epsilon}_s');
p4.Color = c(2,:); % plot meridional strain rate
p5 = plot(rv(1,3:end), fliplr(eps_se), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of elastic deformation to meridional strain rate \dot{\epsilon}_{se}');
p5.Color = c(3,:); % plot \dot{\epsilon}_{se}
p6 = plot(rv(1,2:end), fliplr(gam_s), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of exocytosis to meridional strain rate \gamma_s');
p6.Color = c(4,:);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
set(ax, 'fontsize', 12);
ax.YColor = 'black';
ylabel('$\dot{\epsilon}_s = \dot{\epsilon}_{s e} + \gamma_s$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxis.Visible = 'off'; % turn off this x-axis
% title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(gca, 'YTick', []);
xlim([2 5.6]);
line([3.2 3.3], [0 0])
yyaxis right;
ymin = min([eps_theta eps_thetae gam_theta]); % adjust the y bounds for the strain rate decomposition plot
ymax = max([eps_theta eps_thetae gam_theta]);
p4 = plot(rv(1,3:end), fliplr(eps_theta), '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain rate \dot{\epsilon}_s');
p4.Color = c(5,:); % plot meridional strain rate
p5 = plot(rv(1,3:end), fliplr(eps_thetae), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of elastic deformation to meridional strain rate \dot{\epsilon}_{se}');
p5.Color = c(6,:); % plot \dot{\epsilon}_{se}
p6 = plot(rv(1,2:end), fliplr(gam_theta), '-', 'LineWidth', 2.0, 'DisplayName', 'contribution of exocytosis to meridional strain rate \gamma_s');
p6.Color = c(7,:);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
set(ax, 'fontsize', 12);
ax.YColor = 'black';
set(ax, 'YDir', 'reverse'); % flip the y-axis upside down
ylabel('$\dot{\epsilon}_\theta  = \dot{\epsilon}_{\theta e} + \gamma_\theta$', 'Interpreter', 'latex');
ax = gca;
set(ax, 'fontsize', 12);
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
ax.XRuler.TickLabelGapOffset = -20;
pause(0.5);
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_both_strain_rate_decompositions_hyphoid_a_', num2str(a), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% % plot for grant report
% hold on;
% c = lines(8);
% xlabel("profile");
% % yyaxis left;
% p1 = plot([rv_init(1,:) fliplr(rv_init(1,:))], [rv_init(2,:) -fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
% p1.Color = c(1,:);
% p2 = plot([rv(1,:) fliplr(rv(1,:))], [rv(2,:) -fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
% p2.Color = 'black';
% p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
% p3.Color = c(2,:);
% p4 = plot(rv(1,2:end), velocity / max(velocity), '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
% p4.Color = c(5,:);
% ylabel("profile, secretion rate, velocity");
% ax = gca;
% ax.YColor = 'black';
% % yyaxis right;
% % p5 = plot(rv(1,2:end), strainl, '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain \lambda_s');
% % p5.Color = c(3,:);
% % p6 = plot(rv(1,2:end), strainr, '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain \lambda_\theta');
% % p6.Color = c(4,:);
% % ylabel("strains");
% xlim([0 max(rv(1,:))]);
% legend([p1 p2 p3 p4], 'Location', 'southwest');
% ax = gca;
% ax.YColor = 'black';
% title("Deformation and growth");
% exportgraphics(gcf, ['media/anisotropic_profile_overlaid_grant_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
% close all;
% 
% % plot for grant report - alternate version
% hold on;
% c = lines(8);
% xlabel("profile");
% yyaxis left;
% p1 = plot([rv_init(1,:) fliplr(rv_init(1,:))], [rv_init(2,:) -fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
% p1.Color = c(1,:);
% p2 = plot([rv(1,:) fliplr(rv(1,:))], [rv(2,:) -fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
% p2.Color = 'black';
% ylabel("profile");
% ax = gca;
% ax.YColor = 'black';
% yyaxis right;
% p3 = plot(rv(1,2:end), fliplr(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
% p3.Color = c(2,:);
% p4 = plot(rv(1,2:end), velocity / max(velocity), '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
% p4.Color = c(5,:);
% ylabel("secretion rate, velocity");
% xlim([0 max(rv(1,:))]);
% legend([p1 p2 p3 p4], 'Location', 'southwest');
% ax = gca;
% ax.YColor = 'black';
% title("Deformation and growth");
% exportgraphics(gcf, ['media/anisotropic_profile_overlaid_grant_alt_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
% close all;

etime = cputime-time

end

function sigmaS = sigmaS(strainS, strainR, k, mu)
    sigmaS = mu/2 * (strainR.^-2 - strainS.^-2) + k * (strainR .* strainS - 1);
end

function sigmaTheta = sigmaTheta(strainS, strainR, k, mu)
    sigmaTheta = mu/2 * (strainS.^-2 - strainR.^-2) + k * (strainR .* strainS - 1);
end