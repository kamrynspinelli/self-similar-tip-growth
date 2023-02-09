function [ Tl,Tr,rv,error ] = equilibrium_newton_parabolic_degenerate( rv,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
    global Tol Rtol TolFun TolX Inc 
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt'); 
    error = 10*Tol; 
    inc = Inc;
    
%Find initial guess near the solution
    while error>Tol 
        [tX,X] = ode45(@solver_parabolic_degenerate,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
        error = max(abs(X(end,:)-X0)) % picking out the maximum component
        % makes it terminate faster for small N
        % error = norm(X(end,:)-X0) % using the norm of all the components 
        % better reflects the error computation
        X0=X(end,:);
        clear X;
        rv = reshape(X0',N,2)';  
        rb = LUT*rv';
        D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
        rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
        angle = atan(rb(:,2)./rb(:,1));
    end
  
    X = fsolve(@solver_parabolic_degenerate_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    X0=X;
    X = fsolve(@solver_parabolic_degenerate_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    error = max(abs(X-X0))
    X0=X;
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
    angle = atan(rb(:,2)./rb(:,1));
    
    currentArcs = ParabolicArc.all_degenerate_arcs(rv, K, MU, L0, R0);
    for i = 1:N-1 % there is one fewer patch than there are marker points
        Tl(i) = currentArcs(i+1).sigmaS;
        Tr(i) = currentArcs(i+1).sigmaTheta;
    end
    Tl = Tl'; % turn into column vectors
    Tr = Tr';
end