clear all
time = cputime;
%% parameters
global Tol Rtol TolFun TolX Inc L P
global Fext external_force_mode variance_factor triangle_half_width

% simulation parameters
model = 'spline_local'; % one of 'linear', 'parabolic', 'degenerate', 'spline', 'spline_midpt', 'spline_local'

% dimensional parameters
L = 1;%characteristic length
P = 2;%characteristic pressure
N = 128; %initial number of element

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
% Inc = 0.03;

% external force parameters
external_force_mode = "triangular"; % one of "tip_only", "gaussian", or "triangular"
Fext = 3;
variance_factor = 1/8;
triangle_half_width = 0.1;

%% Initialize cells
R = 1;%the intrinsic radius, if a sphere is to be generated
L0Linear = ones(N,1);
R0Linear = ones(N,1);
L0Parabolic = ones(N,1);
R0Parabolic = ones(N,1);
[ rv, adj ] = generate_sphere(R,N);
% rv(1,:) = 2 * rv(1,:); % elliptic shell % ===== DEBUG =====
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
[L0Spline R0Spline] = spline_intrinsics(rv);
L0SplineLocal = [0 cumsum(L0Spline)'];
R0SplineLocal = rv(2,:);
KSplineLocal = k * ones(N+1, 1);
MUSplineLocal = mu * ones(N+1, 1);
[L0SplineMidpt lR0SplineMidpt rR0SplineMidpt] = spline_intrinsics_strain_midpt(rv);

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

% Solve without growth events
ext_verts = find(sum(adj,1)==1);

if model == "parabolic"
    % [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline(rv,LUT,L0Spline,R0Spline,K,MU,ext_verts,0); % parabolic arcs
elseif model == "spline_midpt"
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_strain_midpt(rv,LUT,L0SplineMidpt,lR0SplineMidpt,rR0SplineMidpt,K,MU,ext_verts,0); % parabolic arcs
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

% Solve with external force
ext_verts = find(sum(adj,1)==1);

if model == "parabolic"
    % [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_parabolic(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,1); % parabolic arcs
    [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,1); % parabolic arcs
elseif model == "linear"
    [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,1); % linear segments
elseif model == "spline"
    [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_cspline(rv,LUT,L0Spline,R0Spline,K,MU,ext_verts,1); % splines
elseif model == "spline_midpt"
    [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_cspline_strain_midpt(rv,LUT,L0SplineMidpt,lR0SplineMidpt,rR0SplineMidpt,K,MU,ext_verts,1); % splines
elseif model == "spline_local"
    [ Tl,Tr,rv,error,frames2,strainFramesl2,strainFramesr2 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,1,L0Parabolic,R0Parabolic,K,MU); % local splines
end

toc;

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

f = [frames1; frames2]; % marker point positions over time
sl = [strainFramesl1; strainFramesl2]; % meridional strains over time
sr = [strainFramesr1; strainFramesr2]; % circumferential strains over time
if model == "linear"
    cmap = jet;
end

% Movie: patches and marker point positions
if external_force_mode == "gaussian"
    v = VideoWriter(['media/timelapse-patches-', model, '-sphere-gaussian-Fext-', num2str(Fext), ...
        '-var-', num2str(variance_factor), '.avi']);
elseif external_force_mode == "triangular"
    v = VideoWriter(['media/timelapse-patches-', model, '-sphere-triangular-Fext-', num2str(Fext), ...
        '-width-', num2str(triangle_half_width), '.avi']);
end
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
    elseif model == "spline" || model == "spline_local"
        slope_cutoff = -1; % get the spline data
        cut = 2;
        slope = (p(2,3) - p(2,1)) / (p(1,3) - p(1,1));
        while slope > slope_cutoff
            cut = cut+1;
            slope = (p(2,cut+1) - p(2,cut-1)) / (p(1,cut+1) - p(1,cut-1));
        end
        [horizA horizB horizC horizD] = find_splines(p(1,1:cut), p(2,1:cut), 0, slope);
        [vertA vertB vertC vertD] = find_splines(-p(2,cut:end), p(1,cut:end), -1/slope, 0);
        hold on; % do the plot
        for i = 1:size(horizA,2)
            spl = @(x) horizA(i) * x.^3 + horizB(i) * x.^2 + horizC(i) * x + horizD(i);
            plot(p(1,i):0.001:p(1,i+1), spl(p(1,i):0.001:p(1,i+1)));
        end
        for i = 1:size(vertA,2)
            spl = @(x) vertA(i) * x.^3 + vertB(i) * x.^2 + vertC(i) * x + vertD(i);
            plot(spl(-p(2,cut+i-1):0.001:-p(2,cut+i)), -(-p(2,cut+i-1):0.001:-p(2,cut+i)));
        end
    end
    daspect([1 1 1]);
    xlim([0 1.5]); ylim([0 1.5]);
    pause(0.2);
    thisFrame = getframe(gcf);
    close all;
    writeVideo(v, thisFrame);
end
close(v);

% Movie: meridional strains
v = VideoWriter(['media/timelapse-strainl-', model, '-sphere-Fext-', num2str(Fext), ...
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
v = VideoWriter(['media/timelapse-strainr-', model, '-sphere-Fext-', num2str(Fext), ...
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
v = VideoWriter(['media/timelapse-Tl-', model, '-sphere-Fext-', num2str(Fext), ...
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
v = VideoWriter(['media/timelapse-Tr-', model, '-sphere-Fext-', num2str(Fext), ...
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