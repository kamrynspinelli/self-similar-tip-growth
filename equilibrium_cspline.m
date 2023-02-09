function [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_newton_spline( rv,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
    global Tol Rtol TolFun TolX Inc 
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt'); 
    error = 10*Tol; 
    inc = Inc;
    % initialize the first frame data
    frames = X0;
    [L R] = spline_intrinsics(rv);
    strainFramesl(1,:) = L ./ L0;
    strainFramesr(1,:) = R ./ R0;
    
%Find initial guess near the solution
    while error>Tol
        [tX,X] = ode45(@solver_cspline,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
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
        sRow = size(strainFramesl, 1) + 1;
        [L R] = spline_intrinsics(rv);
        strainFramesl(sRow,:) = L ./ L0;
        strainFramesr(sRow,:) = R ./ R0;
    end
  
    X = fsolve(@solver_cspline_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    sRow = size(strainFramesl, 1) + 1;
    [L R] = spline_intrinsics(rv);
    strainFramesl(sRow,:) = L ./ L0;
    strainFramesr(sRow,:) = R ./ R0;
    X = fsolve(@solver_cspline_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % parabolic arcs
    error = max(abs(X-X0))
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    sRow = size(strainFramesl, 1) + 1;
    [L R] = spline_intrinsics(rv);
    strainFramesl(sRow,:) = L ./ L0;
    strainFramesr(sRow,:) = R ./ R0;
    rv = reshape(X0',N,2)';  
%     rb = LUT*rv';
%     D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
%     rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
%     angle = atan(rb(:,2)./rb(:,1));
    Tl = sigmaS(strainFramesl(end,:), strainFramesr(end,:), K', MU');
    Tr = sigmaTheta(strainFramesl(end,:), strainFramesr(end,:), K', MU');
    Tl = Tl'; % turn into column vectors
    Tr = Tr';
end

function sigmaS = sigmaS(strainS, strainR, k, mu)
    sigmaS = mu/2 .* (strainR.^-2 - strainS.^-2) + k .* (strainR .* strainS - 1);
end

function sigmaTheta = sigmaTheta(strainS, strainR, k, mu)
    sigmaTheta = mu/2 .* (strainS.^-2 - strainR.^-2) + k .* (strainR .* strainS - 1);
end