function [] = growth_anisotropic_ellipse(a, b, stiffness)
% clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% a = 1;
% b = 3;

% simulation parameters
model = 'parabolic'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

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
Inc = 0.25;
% Inc = 0.05; % slow movie speed


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
gam = zeros(1, size(L0Parabolic, 2) + 1);
gam_s = zeros(1, size(L0Parabolic, 2) + 1);
gam_theta = zeros(1, size(L0Parabolic, 2) + 1);
s(1) = 0;
gam(1) = 1;
% gam_s(1) = 1;
% gam_theta(1) = 1;
gam_s(1) = strainl(end) - 1;
gam_theta(1) = strainl(end) - 1;
for i = 1:size(L0Parabolic, 1)-1
    s(i+1) = s(i) + L0Parabolic(N-i+1);
%     gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
%         * (strainl(N-i+1) / strainr(N-i+1)) / rv_init(2,N-i+1) ...
%         - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
%         * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
%     gam_theta(i+1) = (((rv(2,N-i+1) - rv(2,N-i+2)) / currentArclength(N-i+1)) ...
%         * (strainl(N-i+1) / rv(2,N-i+1)) ...
%         - ((strainr(N-i) - strainr(N-i+1)) / L0Parabolic(N-i+1)) / strainr(N-i+1)) ...
%         * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
    gam_theta(i+1) = (rv(2,N-i+1) - rv(2,N-i+2)) / L0Parabolic(N-i+1) / rv(2,N-i+1) ...
        * gam_s(1:i) * L0Parabolic(N:-1:N-i+1);
    gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
    gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
end
velocity = fliplr(cumsum(s .* gam_s));

% plot \gamma_s and \gamma_\theta
plot(s, gam_s, s, gam_theta, 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Anisotropic growth profile: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min([gam_s gam_theta]) max([gam_s gam_theta])]);
xlabel("s");
ylabel("\gamma_s(s), \gamma_\theta(s)");
exportgraphics(gcf, ['media/anisotropic_growth_gamma_s_theta_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
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
exportgraphics(gcf, ['media/anisotropic_growth_gamma_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot strains
plot(s, fliplr(strainl), s, fliplr(strainr), 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Strain distributions: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min([strainl strainr]) max([strainl strainr])]);
xlabel("s");
ylabel("\lambda_s(s), \lambda_\theta(s)");
exportgraphics(gcf, ['media/anisotropic_strain_distribution_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot monotonicity condition
r = fliplr(rv_init(2,:));
r_p = diff(r) ./ fliplr(L0Parabolic');
r_pp = diff(r_p) ./ fliplr(L0Parabolic(1:end-1)');
strain_tilde = fliplr((strainl - 1) ./ (strainr - 1));
strain_theta_star = fliplr(strainr - 1);
strain_theta_star_p = diff(strain_theta_star) ./ fliplr(L0Parabolic(1:end-1)');
Y = r_pp + (strain_tilde(1:end-1) - 1) .* r_p(1:end-1).^2 ./ r(1:end-2) - strain_theta_star_p ./ strain_theta_star(1:end-1) .* r_p(1:end-1);
hold on;
plot(s(2:end), Y, 'LineWidth', 2.0);
title("Monotonicity of \gamma")
xlabel("s")
ylabel("(r^0)'' + (\lambda - 1) ((r^0)')^2 / r^0 - \lambda_\theta' / \lambda_\theta^* (r^0)'");
p = plot(s, zeros(size(s)), '--', 'LineWidth', 2.0);
p.Color = 'black';
exportgraphics(gcf, ['media/anisotropic_monotonicity_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
close all;

% plot all quantities wrt z
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
ax.YColor = 'black';
xlim([0 3.2]);
yyaxis right;
p3 = plot(rv(1,2:end), fliplr(gam) / max(gam), '-', 'LineWidth', 2.0, 'DisplayName', 'non-dimensional secretion rate \gamma');
p3.Color = c(2,:); % plot \gamma, normalized by its maximum
p4 = plot(rv(1,2:end), velocity, '-', 'LineWidth', 2.0, 'DisplayName', 'rescaled tangential |v|');
p4.Color = c(5,:); % plot velocity
ylim([0 1.1*max([gam/max(gam) velocity])]); % adjust the y-bounds
ax.YColor = 'black';
ylabel("secretion rate, velocity");
ax = gca;
ax.XAxis.Visible = 'off'; % turn off this x-axis
title("Deformation, stretch ratios, and growth");
nexttile; % lower panel
hold on;
yyaxis left;
p1 = plot([fliplr(rv_init(1,:))], [-fliplr(rv_init(2,:))], '--', 'LineWidth', 2.0, 'DisplayName', 'profile without turgor pressure');
p1.Color = c(1,:); % plot lower unturgored profile
p2 = plot([fliplr(rv(1,:))], [-fliplr(rv(2,:))], '-', 'LineWidth', 2.0, 'DisplayName', 'profile with turgor pressure');
p2.Color = 'black'; % plot lower turgored profile
ylim([-1.2 0]);
ax = gca;
ax.YColor = 'black';
xlim([0 3.2]);
line([3.2 3.3], [0 0])
yyaxis right;
p5 = plot(rv(1,2:end), strainl, '-', 'LineWidth', 2.0, 'DisplayName', 'meridional strain \lambda_s');
p5.Color = c(3,:); % plot meridional strain
p6 = plot(rv(1,2:end), strainr, '-', 'LineWidth', 2.0, 'DisplayName', 'circumferential strain \lambda_\theta');
p6.Color = c(4,:); % plot circumferential strain
ymin = min([strainl strainr]); % adjust the y bounds for the strain plot
ymax = max([strainl strainr]);
ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ax = gca;
ax.YColor = 'black';
ylabel("stretch ratios");
ax = gca;
ax.XAxisLocation = 'top'; % move the x-axis for the lower plot to its top
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_ellipse_a_', num2str(a), '_b_', num2str(b), '_k_', num2str(k), '_mu_', num2str(mu), '.png']);
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