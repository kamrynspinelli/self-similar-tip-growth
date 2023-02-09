clear all;

% universal parameters
global Tol Rtol TolFun TolX Inc L P
P = 1; % pressure
Rtol = 1e-13; % computational parameters
Tol =  1e-3;
TolFun = 1e-16;
TolX = 1e-16;
Inc = 0.05;
k = 1; % material property
mu = 1;
ext_force_status = 0;
options = odeset('RelTol',Rtol);

steady_r = (1 + sqrt(17)) / 4; % steady-state radius

% discretization bounds
pow_min = 4;
pow_max = 8;

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
    for i=N:-1:1
        initialArcs = ParabolicArc.all_arcs(rv, ones(1,N), ones(1,N), ones(1,N), ones(1,N));
        L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
        R0Parabolic(i) = initialArcs(i+1).vert(2);
    end
    KParabolic = k * ones(N, 1);
    MUParabolic = mu * ones(N, 1);
    [L0Spline R0Spline] = spline_intrinsics(rv);
    L0SplineLocal = [0 cumsum(L0Spline)'];
    R0SplineLocal = rv(2,:);
    L0SplineLocalInit = L0SplineLocal;
    R0SplineLocalInit = R0SplineLocal;
    KSplineLocal = k * ones(N+1, 1);
    MUSplineLocal = mu * ones(N+1, 1);
    
    % save the true steady state
    rv_steady = steady_r * rv;
    
    % find the initial guess for the solver
    error = 10*Tol; 
    X0 = reshape(rv',2*(N+1),1)';
    while error>Tol 
%         [tX,X] = ode45(@solver_cspline_local,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts,ext_force_status);
        [tX,X] = ode45(@solver_parabolic_imperative,[0 Inc],X0,options,LUT,L0Parabolic,R0Parabolic,KParabolic,MUParabolic,ext_verts,ext_force_status); % get initial guess using the linear solver
        error = max(abs(X(end,:)-X0)) % picking out the maximum component
        % makes it terminate faster for small N
        % error = norm(X(end,:)-X0) % using the norm of all the components 
        % better reflects the error computation
        X0=X(end,:);
        clear X;
        rv = reshape(X0',N+1,2)'; 
    end
    
    % run the solver and find the steady state
    rv = run_nonlinear_solver(rv, LUT, L0SplineLocal, R0SplineLocal, KSplineLocal, MUSplineLocal, ext_verts);
    rvs{exp} = rv;
    mat = rv - rv_steady;
    err1(exp) = norm(mat(:), 1);
    err2(exp) = norm(mat(:), 2);
    errInf(exp) = norm(mat(:), Inf);
end

loglog(2.^pow, err1, 2.^pow, err2, 2.^pow, errInf);
xlabel('N');
ylabel('positional error');
legend(['1-norm, slope = ', num2str(mean(diff(log(err1)) ./ diff(log(2.^pow))))], ...
    ['2-norm, slope = ', num2str(mean(diff(log(err2)) ./ diff(log(2.^pow))))], ...
    ['\infty-norm, slope = ', num2str(mean(diff(log(errInf)) ./ diff(log(2.^pow))))]);
exportgraphics(gcf, 'media/positional-error-analysis-sphere.png')

function res = run_nonlinear_solver(rv, LUT, L0, R0, K, MU, ext_verts)
    global TolFun TolX
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options2 = optimoptions('fsolve','TolFun',TolFun * 4,'TolX',TolX * 4,'Algorithm','levenberg-marquardt');
%     options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX);
%     [X fval outflag] = fsolve(@solver_cspline_local_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,0);
    X = newt_vec(@solver_cspline_local_fast, X0, TolFun, TolX * 4, LUT,L0,R0,K,MU,ext_verts,0);
    res = reshape(X',N,2)'; 
end