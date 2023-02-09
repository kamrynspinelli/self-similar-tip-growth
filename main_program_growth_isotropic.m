clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor

% simulation parameters
model = 'linear'; % one of 'linear', 'parabolic', 'degenerate'

% dimensional parameters
L = 1;%characteristic length
P = 1;%characteristic pressure
N = 16; %initial number of element

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

% external force parameters
external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
Fext = 0;
variance_factor = 1/4;

% growth parameters
growth_tstep = 0.05;
tip_growth_rate = 1;

%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated
L0Linear = ones(N,1);
R0Linear = ones(N,1);
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);
% spherical configuration
[ rv, adj ] = generate_sphere(R,N);
% tubular configuration
% [ rv, adj ] = generate_sphere(R,N/2);
% rv(1,:) = rv(1,:) + 1; % translate the quarter-circle to the right
% rv(:,N/2+1:N+1) = rv(:,1:N/2+1); % make space for the beginning of the array
% for i = 0:N/2-1
%     rv(:,i+1) = [2*i/N; R];
% end
% [ blah, adj ] = generate_sphere(R, N); % set up the connectivity data
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
for i=N:-1:1
   index = find(LUT(i,:));
   L0Linear(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2); % linear segments
   R0Linear(i) = (rv(2,index(1))+rv(2,index(2)))/2;
   L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
   R0Parabolic(i) = initialArcs(i+1).vert(2);
end

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

for i=1:N % save the cos(alpha_0) data for later
    cos_angle(i) = (rv(2,i+1) - rv(2,i)) / (sqrt((rv(2,i+1) - rv(2,i))^2 + (rv(1,i+1) - rv(1,i))^2));
    sin_angle(i) = (rv(1,i+1) - rv(1,i)) / (sqrt((rv(2,i+1) - rv(2,i))^2 + (rv(1,i+1) - rv(1,i))^2));
end

% Solve without growth events
ext_verts = find(sum(adj,1)==1);

if model == "parabolic"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
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

% Solve with external force
if model == "parabolic"
    [ Tl,Tr,rv,error,framestmp,strainFramesltmp,strainFramesrtmp ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,1); % parabolic arcs
    frames = [frames; framestmp];
    strainFramesl = [strainFramesl; strainFramesltmp];
    strainFramesr = [strainFramesr; strainFramesrtmp];
elseif model == "linear"
    [ Tl,Tr,rv,error,framestmp,strainFramesltmp,strainFramesrtmp ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,1); % linear segments
    frames = [frames; framestmp];
    strainFramesl = [strainFramesl; strainFramesltmp];
    strainFramesr = [strainFramesr; strainFramesrtmp];
end

% Solve with growth
L0ParabolicInitial = L0Parabolic;
R0ParabolicInitial = R0Parabolic;
L0LinearInitial = L0Linear;
R0LinearInitial = R0Linear;
for growth_inc = 1:30
    if model == "parabolic"
        L0Parabolic = L0ParabolicInitial .* exp(i/5 * tip_growth_rate * -cos_angle)';
        R0Parabolic = R0ParabolicInitial .* exp(i/5 * tip_growth_rate * -cos_angle)';
        [ Tl,Tr,rv,error,framestmp,strainFramesltmp,strainFramesrtmp ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,1); % parabolic arcs
        frames = [frames; framestmp];
        strainFramesl = [strainFramesl; strainFramesltmp];
        strainFramesr = [strainFramesr; strainFramesrtmp];
    elseif model == "linear"
        a = 0;
        for i=1:N % find the cos(alpha_0) data corresponding to the current L0/R0 configuration
            cos_angle(N-i+1) = 2 * (a - R0Linear(N-i+1)) / L0Linear(N-i+1);
            a = a - L0Linear(N-i+1) * cos_angle(N-i+1);
        end
        if max(abs(cos_angle)) > 1
            1;
        end
%         L0Linear = L0LinearInitial .* exp(i/5 * tip_growth_rate * cos_angle)';
%         R0Linear = R0LinearInitial .* exp(i/5 * tip_growth_rate * cos_angle)';
        L0Linear = L0Linear * tip_growth_rate * growth_tstep .* -cos_angle' + L0Linear;
        R0Linear = R0Linear * tip_growth_rate * growth_tstep .* -cos_angle' + R0Linear;
        [ Tl,Tr,rv,error,framestmp,strainFramesltmp,strainFramesrtmp ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,1); % linear segments
        frames = [frames; framestmp];
        strainFramesl = [strainFramesl; strainFramesltmp];
        strainFramesr = [strainFramesr; strainFramesrtmp];
    end
end

Eg(1) = error;
%record initial condition
[cl,cr] = compute_curvatures(rv);
cl_t(3).dat = cl;
cr_t(3).dat = cr;
rv_t(3).dat = rv;
LUT_t(3).dat = LUT;
adj_t(3).dat = adj;
L0_t(3).dat = L0Parabolic;
R0_t(3).dat = R0Parabolic;
Tl_t(3).dat = Tl;
Tr_t(3).dat = Tr;

%%

% [Tmax Tmin]=plot_fig_t(rv_t,LUT_t,Tl_t,Tr_t,R0Parabolic,K,'tension');
% plot_fig_t(rv_t,LUT_t,cl_t,cr_t,R0Parabolic,K,'curvature');

f = frames; % marker point positions over time
sl = strainFramesl; % meridional strains over time
sr = strainFramesr; % circumferential strains over time
if model == "linear"
    cmap = jet;
end

% Movie: patches and marker point positions
v = VideoWriter(['media/timelapse-growth-isotropic-patches-', model, '-sphere-Fext-', num2str(Fext), ...
    '-var-', num2str(variance_factor), '.avi']);
open(v);
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    if model == "parabolic"
        ParabolicArc.plot_arcs(p);
    elseif model == "linear"
        hold on;
        for i=1:N
            plot(p(1,i:i+1), p(2,i:i+1), 'LineWidth', 2.0);
        end
    end
    daspect([1 1 1]);
    xlim([0 4.5]); ylim([0 4.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

% ========= DEBUG =========
% Rsteady = 1.0025;
% ParabolicArc.plot_arcs(p);
% X = rv(1,end)-2*Rsteady:0.01:rv(1,end);
% Y = sqrt(Rsteady^2 - (X - p(1,end) + Rsteady).^2);
% hold on;
% plot(X, Y);
% ========= END DEBUG =========


% Movie: meridional strains
v = VideoWriter(['media/timelapse-growth-isotropic-strainl-', model, '-sphere-Fext-', num2str(Fext), ...
    '-var-', num2str(variance_factor), '.avi']);
open(v);
heatMin = min(min(sl));
heatMax = max(max(sl));
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    if model == "parabolic"
        ParabolicArc.plot_arcs(p, sl(frameNum,:), heatMin, heatMax);
    elseif model == "linear"
        hold on;
        for i=1:N
            pl = plot(p(1,i:i+1), p(2,i:i+1), 'LineWidth', 2.0);
            if heatMin == heatMax
                colorVal = 1;
            else
                colorVal = min(max(round((sl(frameNum,i) - heatMin) / (heatMax - heatMin) * 256), 1), 256);
            end
            pl.Color = cmap(colorVal, :);
        end
    end
    daspect([1 1 1]);
    xlim([0 1.5]); ylim([0 1.5]);
    colormap jet; colorbar; caxis([heatMin heatMax]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

% Movie: circumferential strains
v = VideoWriter(['media/timelapse-growth-isotropic-strainr-', model, '-sphere-Fext-', num2str(Fext), ...
    '-var-', num2str(variance_factor), '.avi']);
open(v);
heatMin = min(min(sr));
heatMax = max(max(sr));
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    if model == "parabolic"
        ParabolicArc.plot_arcs(p, sr(frameNum,:), heatMin, heatMax);
    elseif model == "linear"
        hold on;
        for i=1:N
            pl = plot(p(1,i:i+1), p(2,i:i+1), 'LineWidth', 2.0);
            if heatMin == heatMax
                colorVal = 1;
            else
                colorVal = min(max(round((sr(frameNum,i) - heatMin) / (heatMax - heatMin) * 256), 1), 256);
            end
            pl.Color = cmap(colorVal, :);
        end
    end
    daspect([1 1 1]);
    xlim([0 1.5]); ylim([0 1.5]);
    colormap jet; colorbar; caxis([heatMin heatMax]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

% Movie: \sigma_s
v = VideoWriter(['media/timelapse-growth-isotropic-Tl-', model, '-sphere-Fext-', num2str(Fext), ...
    '-var-', num2str(variance_factor), '.avi']);
open(v);
for i = 1:size(sl, 1)
    for j = 1:size(sl,2)
        TlFrames(i,j) = sigmaS(sl(i,j), sr(i,j), k, mu);
    end
end
heatMin = min(min(TlFrames));
heatMax = max(max(TlFrames));
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    if model == "parabolic"
        ParabolicArc.plot_arcs(p, TlFrames(frameNum,:), heatMin, heatMax);
    elseif model == "linear"
        hold on;
        for i=1:N
            pl = plot(p(1,i:i+1), p(2,i:i+1), 'LineWidth', 2.0);
            if heatMin == heatMax
                colorVal = 1;
            else
                colorVal = min(max(round((TlFrames(frameNum,i) - heatMin) / (heatMax - heatMin) * 256), 1), 256);
            end
            pl.Color = cmap(colorVal, :);
        end
    end
    daspect([1 1 1]);
    xlim([0 1.5]); ylim([0 1.5]);
    colormap jet; colorbar; caxis([heatMin heatMax]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

% Movie: \sigma_\theta
v = VideoWriter(['media/timelapse-growth-isotropic-Tr-', model, '-sphere-Fext-', num2str(Fext), ...
    '-var-', num2str(variance_factor), '.avi']);
open(v);
for i = 1:size(sl, 1)
    for j = 1:size(sl,2)
        TrFrames(i,j) = sigmaTheta(sl(i,j), sr(i,j), k, mu);
    end
end
heatMin = min(min(TrFrames));
heatMax = max(max(TrFrames));
for frameNum=1:size(f,1)
    p = reshape(f(frameNum,:), N+1, 2)';
    if model == "parabolic"
        ParabolicArc.plot_arcs(p, TrFrames(frameNum,:), heatMin, heatMax);
    elseif model == "linear"
        hold on;
        for i=1:N
            pl = plot(p(1,i:i+1), p(2,i:i+1), 'LineWidth', 2.0);
            if heatMin == heatMax
                colorVal = 1;
            else
                colorVal = min(max(round((TrFrames(frameNum,i) - heatMin) / (heatMax - heatMin) * 256), 1), 256);
            end
            pl.Color = cmap(colorVal, :);
        end
    end
    daspect([1 1 1]);
    xlim([0 1.5]); ylim([0 1.5]);
    colormap jet; colorbar; caxis([heatMin heatMax]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

close all;
etime = cputime-time

function sigmaS = sigmaS(strainS, strainR, k, mu)
    sigmaS = mu/2 * (strainR^-2 - strainS^-2) + k * (strainR * strainS - 1);
end

function sigmaTheta = sigmaTheta(strainS, strainR, k, mu)
    sigmaTheta = mu/2 * (strainS^-2 - strainR^-2) + k * (strainR * strainS - 1);
end

function dGamma = dGamma(s0, arcs, initial_arcs, L0, R0, cosAlpha0, sinAlpha0) % (dGamma/ds^0) / Gamma
%     usage: 
%     arcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), L0Parabolic, R0Parabolic);
%     arcs = arcs(2:end-1); % get current arcs
%     dGamma(2, arcs, initialArcs, L0Parabolic, R0Parabolic, cos_angle) % value of (dGamma/ds^0) / Gamma
%     [s,y] = ode45(@(s,y) y * dGamma(s, arcs, initialArcs, L0Parabolic, R0Parabolic, cos_angle), [0 2], 1) % solve the ODE
    patch = size(arcs, 2); % the last patch
    % arclengthSum = initial_arcs(patch).arclength;
    arclengthSum = L0(patch);
    while arclengthSum < s0 % find which patch this s-coordinate corresponds to
        patch = patch - 1;
        % arclengthSum = arclengthSum + initial_arcs(patch).arclength;
        arclengthSum = arclengthSum + L0(patch);
    end
    dGamma = -cosAlpha0(patch) / R0(patch);
end