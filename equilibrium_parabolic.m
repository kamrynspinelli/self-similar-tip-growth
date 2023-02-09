function [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_newton_parabolic( rv,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
    global Tol Rtol TolFun TolX Inc 
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt'); 
    error = 10*Tol; 
    inc = Inc;
    % initialize the first frame data
    frames = X0;
    currentArcs = ParabolicArc.all_arcs(rv, K, MU, L0, R0);
    for i=1:N-1
        strainFramesl(1,i) = currentArcs(i+1).arclength / L0(i);
        strainFramesr(1,i) = currentArcs(i+1).vert(2) / R0(i);
    end
    
%Find initial guess near the solution
    while error>Tol 
        [tX,X] = ode45(@solver_parabolic,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
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
        frames(end+1, :) = X0;
        currentArcs = ParabolicArc.all_arcs(rv, K, MU, L0, R0);
        sRow = size(strainFramesl, 1) + 1;
        for i=1:N-1
            strainFramesl(sRow,i) = currentArcs(i+1).arclength / L0(i);
            strainFramesr(sRow,i) = currentArcs(i+1).vert(2) / R0(i);
        end
    end
  
    X = fsolve(@solver_parabolic_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    currentArcs = ParabolicArc.all_arcs(rv, K, MU, L0, R0);
    sRow = size(strainFramesl, 1) + 1;
    for i=1:N-1
        strainFramesl(sRow,i) = currentArcs(i+1).arclength / L0(i);
        strainFramesr(sRow,i) = currentArcs(i+1).vert(2) / R0(i);
    end
    X = fsolve(@solver_parabolic_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    error = max(abs(X-X0))
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    currentArcs = ParabolicArc.all_arcs(rv, K, MU, L0, R0);
    sRow = size(strainFramesl, 1) + 1;
    for i=1:N-1
        strainFramesl(sRow,i) = currentArcs(i+1).arclength / L0(i);
        strainFramesr(sRow,i) = currentArcs(i+1).vert(2) / R0(i);
    end
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
    angle = atan(rb(:,2)./rb(:,1));
    
    currentArcs = ParabolicArc.all_arcs(rv, K, MU, L0, R0);
    for i = 1:N-1 % there is one fewer patch than there are marker points
        Tl(i) = currentArcs(i+1).sigmaS;
        Tr(i) = currentArcs(i+1).sigmaTheta;
    end
    Tl = Tl'; % turn into column vectors
    Tr = Tr';
end