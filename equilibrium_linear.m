function [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_newton_linear( rv,LUT,L0,R0,K,MU,ext_verts,ext_force_status)
    global  Tol Rtol TolFun TolX Inc 
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX); 
    error = 10*Tol; 
    inc = Inc;
    % initialize the first frame data
    frames = X0;
    rb = LUT*rv';
    strainFramesl = (sqrt(sum(rb.^2,2)) ./ L0)'; % vector of edge lengths / intrinsic lengths
    strainFramesr = (0.5*(rv(2,1:end-1)+rv(2,2:end))' ./ R0)'; % vector of patch radii / intrinsic radii
    
%Find initial guess near the solution
    while error>Tol 
        [tX,X] = ode45(@solver,[0 inc],X0,options,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % linear segments
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
        strainFramesl(end+1,:) = (D ./ L0)'; % vector of edge lengths / intrinsic lengths
        strainFramesr(end+1,:) = (rm ./ R0)'; % vector of patch radii / intrinsic radii
    end

	X = fsolve(@solver_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % linear segments
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    strainFramesl(end+1,:) = (sqrt(sum(rb.^2,2)) ./ L0)'; % vector of edge lengths / intrinsic lengths
    strainFramesr(end+1,:) = (0.5*(rv(2,1:end-1)+rv(2,2:end))' ./ R0)'; % vector of patch radii / intrinsic radii
	X = fsolve(@solver_fast,X0,options2,LUT,L0,R0,K,MU,ext_verts,ext_force_status); % linear segments
    error = max(abs(X-X0))
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    strainFramesl(end+1,:) = (sqrt(sum(rb.^2,2)) ./ L0)'; % vector of edge lengths / intrinsic lengths
    strainFramesr(end+1,:) = (0.5*(rv(2,1:end-1)+rv(2,2:end))' ./ R0)'; % vector of patch radii / intrinsic radii
    rv = reshape(X0',N,2)';  
    rb = LUT*rv';
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
    angle = atan(rb(:,2)./rb(:,1));

    Tl = K.*((D./L0.*rm./R0)-1)+0.5*MU.*((R0.^2./rm.^2)-L0.^2./D.^2); % linear segments
    Tr = K.*((D./L0.*rm./R0)-1)+0.5*MU.*(L0.^2./D.^2-(R0.^2./rm.^2));
end