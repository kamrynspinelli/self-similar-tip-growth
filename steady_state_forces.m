%% Comparison between the forces with different solvers in the true steady state
function [fLinear, fParabolic, fParabolicDegenerate, fSpline, fSplineLocal] =  steady_state_forces(patches, k, mu)
    % patches: number of patches
    % k: bulk modulus
    % mu: shear modulus
    %% parameters
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor

    % dimensional parameters
    L = 1; % characteristic length
    P = 1; % characteristic pressure
    N = patches + 1; % number of marker points

    % uniform material property
    K = k*ones(N-1,1);
    MU = mu*ones(N-1,1);

    % computational parameters
    Rtol = 1e-13;
    Tol =  1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    Inc =  1;

    % set up the steady-state point configuration
    rSteady = 1 / 4 * (1 + sqrt(17));
    rv(1,:) = rSteady * cos(linspace(pi/2,0,N));
    rv(2,:) = rSteady * sin(linspace(pi/2,0,N)); 
    X0 = reshape(rv',2*N,1)';

    % characteristic arclengths and radii
    R = 1;%the intrinsic radius, if a sphere is to be generated
    L0Parabolic = ones(N-1,1);
    L0Linear = ones(N-1,1);
    R0Parabolic = ones(N-1,1);
    R0Linear = ones(N-1,1);
    [ rv, adj ] = generate_sphere(R,N-1);
    initialArcs = ParabolicArc.all_arcs(rv, ones(1,N-1), ones(1,N-1), ones(1,N-1), ones(1,N-1));
    LUT = zeros(N-1,size(adj,1));
    ii = 1;
    for n = 1:N-1
        nv = find(adj(n,:));
        nv = nv(nv>n);
        for nn = nv
            LUT(ii,n) = 1; LUT(ii,nn) = -1;
            ii = ii + 1;
        end
    end
    for i=N-1:-1:1
       index = find(LUT(i,:));
       L0Linear(i) = sqrt((rv(1,index(1))-rv(1,index(2))).^2+(rv(2,index(1))-rv(2,index(2))).^2); % linear segments
       R0Linear(i) = (rv(2,index(1))+rv(2,index(2)))/2;
       L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
       R0Parabolic(i) = initialArcs(i+1).vert(2);
    end
    [L0Spline R0Spline] = spline_intrinsics(rv);
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal = rv(2,:);
    KSplineLocal = k * ones(N, 1);
    MUSplineLocal = mu * ones(N, 1);
    [L0SplineMidpt lR0SplineMidpt rR0SplineMidpt] = spline_intrinsics_strain_midpt(rv);
    ext_verts = find(sum(adj,1)==1);

    %% Force computation
    rSteady = 1 / 4 * (1 + sqrt(17));
    fLinear = reshape((solver_fast(X0,LUT,L0Linear,R0Linear,K,MU,ext_verts,0))', N, 2)';
    fParabolic = reshape((solver_parabolic_imperative(0,X0,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0))', N, 2)';
    % fParabolic = reshape((solver_parabolic_fast(X0,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0))', N, 2)';
    % fParabolic = reshape((solver_parabolic_analytic(0,X0,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0))', N, 2)';
    fParabolicDegenerate = reshape((solver_parabolic_degenerate_fast(X0,LUT,L0Linear,R0Linear,K,MU,ext_verts,0))', N, 2)';
    fSpline = reshape((solver_cspline(1,X0,LUT,L0Spline,R0Spline,K,MU,ext_verts,0))', N, 2)';
    fSplineMidpt = reshape((solver_cspline_strain_midpt(1,X0,LUT,L0SplineMidpt,lR0SplineMidpt,rR0SplineMidpt,K,MU,ext_verts,0))', N, 2)';
    fSplineLocal = reshape((solver_cspline_local(1,X0,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0))', N, 2)';
    % fSplineLocal = (rSteady * pi/2/patches) * reshape((solver_cspline_local(1,X0,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0))', N, 2)';
    % need factor of rSteady * pi/2/patches to convert from pointwise force
    % to integral force
    fLinearNorm = vecnorm(fLinear);
    fParabolicNorm = vecnorm(fParabolic);
    fParabolicDegenerateNorm = vecnorm(fParabolicDegenerate);
    fSplineNorm = vecnorm(fSpline);
    fSplineMidptNorm = vecnorm(fSplineMidpt);
    fSplineLocalNorm = vecnorm(fSplineLocal);
    maxFNorm = max([vecnorm(fLinear) vecnorm(fParabolic) vecnorm(fParabolicDegenerate) vecnorm(fSpline)]);

    %% Plotting
    cmap = parula; % select the colormap for the force magnitudes
    rv(1,:) = rSteady * cos(linspace(pi/2,0,N));
    rv(2,:) = rSteady * sin(linspace(pi/2,0,N)); 
    
    % force magnitude plot: parabolic arcs
    hold on;
    % first patch
    p = plot([rv(1,1) (rv(1,1) + rv(1,2))/2], [rv(2,1) (rv(2,1) + rv(2,2))/2]);
    colorVal = max(round(fParabolicNorm(1) / max(fParabolicNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    for i = 2:N-1 % interior patches
        p = plot([(rv(1,i-1) + rv(1,i))/2   rv(1,i)   (rv(1,i) + rv(1,i+1))/2], [(rv(2,i-1) + rv(2,i))/2   rv(2,i)   (rv(2,i) + rv(2,i+1))/2]);
        colorVal = max(round(fParabolicNorm(i) / max(fParabolicNorm) * 256), 1);
        p.Color = cmap(colorVal, :);
        p.LineWidth = 5;
    end
    % last patch
    p = plot([(rv(1,N-1) + rv(1,N))/2 rv(1,N)], [(rv(2,N-1) + rv(2,N))/2 rv(2,N)]);
    colorVal = max(round(fParabolicNorm(N) / max(fParabolicNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    title(['Force distribution for parabolic arcs: K = ', num2str(k), ', mu = ', num2str(mu)]);
    colormap parula; colorbar; caxis([0 max(fParabolicNorm)]);
    daspect([1 1 1]); pbaspect([1 1 1]);
    %saveas(p, ['media/parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png']);
    exportgraphics(gcf, ['media/parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
    
    % force vector plot: parabolic arcs
    hold on;
    plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);
    quiver(rv(1,:), rv(2,:), fParabolic(1,:), fParabolic(2,:), 'LineWidth', 2.0);
    daspect([1 1 1]);
    exportgraphics(gcf, ['media/force_vectors_parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;

    % force magnitude plot: linear segments
    hold on;
    % first patch
    p = plot([rv(1,1) (rv(1,1) + rv(1,2))/2], [rv(2,1) (rv(2,1) + rv(2,2))/2]);
    colorVal = max(round(fLinearNorm(1) / max(fLinearNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    for i = 2:N-1 % interior patches
        p = plot([(rv(1,i-1) + rv(1,i))/2   rv(1,i)   (rv(1,i) + rv(1,i+1))/2], [(rv(2,i-1) + rv(2,i))/2   rv(2,i)   (rv(2,i) + rv(2,i+1))/2]);
        colorVal = max(round(fLinearNorm(i) / max(fLinearNorm) * 256), 1);
        p.Color = cmap(colorVal, :);
        p.LineWidth = 5;
    end
    % last patch
    p = plot([(rv(1,N-1) + rv(1,N))/2 rv(1,N)], [(rv(2,N-1) + rv(2,N))/2 rv(2,N)]);
    colorVal = max(round(fLinearNorm(N) / max(fLinearNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    title(['Force distribution for linear segments: K = ', num2str(k), ', mu = ', num2str(mu)]);
    colormap parula; colorbar; caxis([0 max(fLinearNorm)]);
    daspect([1 1 1]); pbaspect([1 1 1]);
    %saveas(p, ['media/linear_k_', num2str(k), '_mu_', num2str(mu), '.png']);
    exportgraphics(gcf, ['media/linear_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
    
    % force vector plot: linear segments
    hold on;
    plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);
    quiver(rv(1,:), rv(2,:), fLinear(1,:), fLinear(2,:), 'LineWidth', 2.0);
    daspect([1 1 1]);
    exportgraphics(gcf, ['media/force_vectors_linear_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;

    % force magnitude plot: degenerate parabolic arcs
    hold on;
    % first patch
    p = plot([rv(1,1) (rv(1,1) + rv(1,2))/2], [rv(2,1) (rv(2,1) + rv(2,2))/2]);
    colorVal = max(round(fParabolicDegenerateNorm(1) / max(fParabolicDegenerateNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    for i = 2:N-1 % interior patches
        p = plot([(rv(1,i-1) + rv(1,i))/2   rv(1,i)   (rv(1,i) + rv(1,i+1))/2], [(rv(2,i-1) + rv(2,i))/2   rv(2,i)   (rv(2,i) + rv(2,i+1))/2]);
        colorVal = max(round(fParabolicDegenerateNorm(i) / max(fParabolicDegenerateNorm) * 256), 1);
        p.Color = cmap(colorVal, :);
        p.LineWidth = 5;
    end
    % last patch
    p = plot([(rv(1,N-1) + rv(1,N))/2 rv(1,N)], [(rv(2,N-1) + rv(2,N))/2 rv(2,N)]);
    colorVal = max(round(fParabolicDegenerateNorm(N) / max(fParabolicDegenerateNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    title(['Force distribution for degenerate parabolic arcs: K = ', num2str(k), ', mu = ', num2str(mu)]);
    colormap parula; colorbar; caxis([0 max(fParabolicDegenerateNorm)]);
    daspect([1 1 1]); pbaspect([1 1 1]);
    %saveas(p, ['media/degenerate_parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png']);
    exportgraphics(gcf, ['media/degenerate_parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
    
    % force vector plot: degenerate parabolic arcs
    hold on;
    plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);
    quiver(rv(1,:), rv(2,:), fParabolicDegenerate(1,:), fParabolicDegenerate(2,:), 'LineWidth', 2.0);
    daspect([1 1 1]);
    exportgraphics(gcf, ['media/force_vectors_degenerate_parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
    
    % force magnitude plot: cubic splines
    hold on;
    % first patch
    p = plot([rv(1,1) (rv(1,1) + rv(1,2))/2], [rv(2,1) (rv(2,1) + rv(2,2))/2]);
    colorVal = max(round(fSplineNorm(1) / max(fSplineNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    for i = 2:N-1 % interior patches
        p = plot([(rv(1,i-1) + rv(1,i))/2   rv(1,i)   (rv(1,i) + rv(1,i+1))/2], [(rv(2,i-1) + rv(2,i))/2   rv(2,i)   (rv(2,i) + rv(2,i+1))/2]);
        colorVal = max(round(fSplineNorm(i) / max(fSplineNorm) * 256), 1);
        p.Color = cmap(colorVal, :);
        p.LineWidth = 5;
    end
    % last patch
    p = plot([(rv(1,N-1) + rv(1,N))/2 rv(1,N)], [(rv(2,N-1) + rv(2,N))/2 rv(2,N)]);
    colorVal = max(round(fSplineNorm(N) / max(fSplineNorm) * 256), 1);
    p.Color = cmap(colorVal, :);
    p.LineWidth = 5;
    title(['Force distribution for cubic splines: K = ', num2str(k), ', mu = ', num2str(mu)]);
    colormap parula; colorbar; caxis([0 max(fSplineNorm)]);
    daspect([1 1 1]); pbaspect([1 1 1]);
    %saveas(p, ['media/degenerate_parabolic_k_', num2str(k), '_mu_', num2str(mu), '.png']);
    exportgraphics(gcf, ['media/spline_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
    
    % force vector plot: cubic splines
    hold on;
    plot(rv(1,:), rv(2,:), 'LineWidth', 2.0);
    quiver(rv(1,:), rv(2,:), fSpline(1,:), fSpline(2,:), 'LineWidth', 2.0);
    daspect([1 1 1]);
    exportgraphics(gcf, ['media/force_vectors_spline_k_', num2str(k), '_mu_', num2str(mu), '.png'],'Resolution',300)
    close all;
end