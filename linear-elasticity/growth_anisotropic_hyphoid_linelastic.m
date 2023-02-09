function [tip_length gam_peak ZN RN Rtip W_int gam L0Spline strainl ksspline] = growth_anisotropic_hyphoid_linelastic(a, poisson, young)
% clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 160; %initial number of element

% Uniform material property
POIS = poisson*ones(N,1);
YOUNG = young*ones(N,1);

%computational parameters
Rtol = 1e-13;
Tol =  1e-3;
TolFun = 1e-16;
TolX = 1e-16;
% Inc =  1; % fast movie speed
% Inc = 0.25;
Inc = 0.05; % slow movie speed
% Inc = 0.01; % molasses movie speed

%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated

[ tmp, adj ] = generate_sphere(R,N/2); % set up the endcap
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
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);
% create ParabolicArc objects for the initial configuration (with dummy tension data)
initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
for i=N:-1:1
   index = find(LUT(i,:));
   L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
   R0Parabolic(i) = initialArcs(i+1).vert(2);
end
[L0Spline R0Spline] = spline_intrinsics(rv);
L0SplineLocal = [0 cumsum(L0Spline)'];
R0SplineLocal = rv(2,:);
POISSplineLocal = poisson * ones(N+1, 1);
YOUNGSplineLocal = young * ones(N+1, 1);

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
L0_t(1).dat = L0SplineLocal;
R0_t(1).dat = R0SplineLocal;
Tl_t(1).dat = Tl;
Tr_t(1).dat = Tr; 

tic;

%% Simulate effects of elastic deformation
ext_verts = find(sum(adj,1)==1);

if model == "spline_local"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local_linelastic(rv,LUT,L0SplineLocal,R0SplineLocal,POISSplineLocal,YOUNGSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,POIS,YOUNG);
end
    
Eg(1) = error;
%record initial condition
[cl,cr] = compute_curvatures(rv);
cl_t(2).dat = cl;
cr_t(2).dat = cr;
rv_t(2).dat = rv;
LUT_t(2).dat = LUT;
adj_t(2).dat = adj;
L0_t(2).dat = L0SplineLocal;
R0_t(2).dat = R0SplineLocal;
Tl_t(2).dat = Tl;
Tr_t(2).dat = Tr;

toc;

%% Induced growth plot
strainl = strainFramesl1(end,:);
strainr = strainFramesr1(end,:);
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
%     gam_theta(i+1) = (rv(2,N-i+1) - rv(2,N-i+2)) / L0Spline(N-i+1) / rv(2,N-i+1) ...
%         * gam_s(1:i) * L0Spline(N:-1:N-i+1);
    gam_theta(i+1) = (rv_init(2,N-i+1) - rv_init(2,N-i+2)) / L0Spline(N-i+1) / rv_init(2,N-i+1) ...
        * gam_s(1:i) * L0Spline(N:-1:N-i+1);
    gam_s(i+1) = gam_theta(i+1) * (strainl(N-i+1) - 1) / (strainr(N-i+1) - 1);
    gam(i+1) = gam_theta(i+1) / (strainr(N-i+1) - 1);
end
velocity = strainl(2:end) .* fliplr([0 cumsum(fliplr(L0Spline(2:end)') .* gam_s(1:end-1))]); % velocity relative to the tip
eps_s = diff(fliplr(velocity)) ./ diff(s); % combined growth rates
eps_theta = 1 ./ fliplr(rv(2,2:end-1)) .* diff(fliplr(rv(2,2:end))) ./ diff(s) .* fliplr(velocity(1:end-1)); % TODO: double-check these indices
[ks ktheta] = compute_curvatures(rv); % estimation by finite difference for each patch
[ksspline kthetaspline] = compute_curvatures_spline(rv); % estimation by cubic splines

s = [0 cumsum(fliplr(L0Spline(2:end)' .* strainl(3:end)))]; % post-deformation arclength

tip_length = strainl(N/2+2:end) * L0Spline(N/2+1:end); % post-deformation arclength
gam_peak = strainl(end-find(gam == max(gam))+2:end) * L0Spline(end-find(gam == max(gam))+2:end);
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

% plot \gamma_s and \gamma_\theta
plot(s, gam_s, s, gam_theta, 'LineWidth', 2.0);
% daspect([1 1 1]);
pbaspect([1 1 1]);
% title(['Anisotropic growth profile: ellipse, a=', num2str(a), ', b = ', num2str(b)]);
xlim([0 max(s)]);
ylim([min([gam_s gam_theta]) max([gam_s gam_theta])]);
xlabel("s");
ylabel("\gamma_s(s), \gamma_\theta(s)");
exportgraphics(gcf, ['media/anisotropic_growth_gamma_s_theta_ellipse_a_', num2str(a), '_pois_', num2str(poisson), '_young_', num2str(young), '.png']);
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
exportgraphics(gcf, ['media/anisotropic_growth_gamma_ellipse_a_', num2str(a), '_pois_', num2str(poisson), '_young_', num2str(young), '.png']);
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
exportgraphics(gcf, ['media/anisotropic_strain_distribution_ellipse_a_', num2str(a), '_pois_', num2str(poisson), '_young_', num2str(young), '.png']);
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
exportgraphics(gcf, ['media/anisotropic_monotonicity_ellipse_a_', num2str(a), '_pois_', num2str(poisson), '_young_', num2str(young), '.png']);
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
xlim([2 rv(1,end)+0.3]);
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
ylim([0 2.9]);
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
xlim([2 rv(1,end)+0.3]);
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
% ylim([ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)]);
ylim([1.004 1.032]);
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
width = 8 * 0.8; height = width * 1.2 / ((max(rv(1,:)) - 2) + 0.3); % width and height for one of the panels
t.InnerPosition = [1 1 width 2*height]; % hack to set the aspect ratio because tiled layout can't do that I guess 
exportgraphics(gcf, ['media/anisotropic_profile_overlaid_hyphoid_a_', num2str(a), '_linelastic_pois_', num2str(poisson), '_young_', num2str(young), '.png']);
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