% universal parameters
global Tol Rtol TolFun TolX Inc L P
P = 1; % pressure
Rtol = 1e-13; % computational parameters
Tol =  1e-3 * 0.7;
TolFun = 1e-16;
TolX = 1e-16;
Inc = 0.05;
k = 1; % material property
mu = 1;

steady_r = (1 + sqrt(17)) / 4; % steady-state radius

% discretization bounds
pow_min = 4;
pow_max = 12;

pow = pow_min:pow_max;
for exp = 1:size(pow,2)
    % set up the unturgid shape
    N = 2^pow(exp); % number of patches
    rv = zeros(2,N+1);
    rv(1,:) = cos(pi/2:-pi/2/N:0);
    rv(2,:) = sin(pi/2:-pi/2/N:0);
    rv(1,1) = 0;
    rv(2,end) = 0;
    
    % set up the distribution of material property
    K = k*ones(N,1);
    MU = mu*ones(N,1);
    
    % make the adjacency table
    LUT = diag(ones(N+1,1)) + diag(-ones(N,1), 1);
    LUT = LUT(1:end-1,:);
    ext_verts = [1 N+1];
    
    % compute the intrinsic data
    [L0Spline R0Spline] = spline_intrinsics(rv);
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal = rv(2,:);
    L0SplineLocalInit = L0SplineLocal;
    R0SplineLocalInit = R0SplineLocal;
    KSplineLocal = k * ones(N+1, 1);
    MUSplineLocal = mu * ones(N+1, 1);
    
    % inflate the shell
    rv = steady_r * rv;
    
    % compute the force residue
    X0 = reshape(rv',2*(N+1),1)';
    F = solver_cspline_local_fast(X0,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0);
    F_r = reshape(F',N+1,2)';
%     hold on; plot(rv(1,:), rv(2,:));
%     quiver(rv(1,:), rv(2,:), F_r(1,:), F_r(2,:));
%     close all;
    res1(exp) = norm(F, 1);
    res2(exp) = norm(F, 2);
    resInf(exp) = norm(F, Inf);
end

loglog(2.^pow, res1, 2.^pow, res2, 2.^pow, resInf);
xlabel('N');
ylabel('force residue');
legend(['1-norm, slope = ', num2str(mean(diff(log(res1)) ./ diff(log(2.^pow))))], ...
    ['2-norm, slope = ', num2str(mean(diff(log(res2)) ./ diff(log(2.^pow))))], ...
    ['\infty-norm, slope = ', num2str(mean(diff(log(resInf)) ./ diff(log(2.^pow))))]);
exportgraphics(gcf, 'media/force-error-analysis-sphere.png')