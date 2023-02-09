function [ Tl,Tr,rv,error,frames,strainFramesl,strainFramesr ] = equilibrium_newton_cspline_local( rv,LUT,L0,R0,POIS,YOUNG,ext_verts,ext_force_status,L0Parabolic,R0Parabolic,POISParabolic,YOUNGParabolic)
    global Tol Rtol TolFun TolX Inc 
    N = size(rv,2);
    X0 = reshape(rv',2*N,1)';
    options = odeset('RelTol',Rtol);
    options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX,'Algorithm','levenberg-marquardt');
    % options2 = optimoptions('fsolve','TolFun',TolFun,'TolX',TolX);
    error = 10*Tol; 
    inc = Inc;
    % initialize the first frame data
    frames = X0;
    [Ls Ltheta] = strains(rv, L0, R0);
    strainFramesl(1,:) = Ls;
    strainFramesr(1,:) = Ltheta;
    
%Find initial guess near the solution
    while error>Tol * 0.8
        % [tX,X] = ode45(@solver_cspline_local_linelastic,[0 inc],X0,options,LUT,L0,R0,POIS,YOUNG,ext_verts,ext_force_status); % parabolic arcs
        [tX,X] = ode45(@solver_parabolic_imperative_linelastic,[0 inc],X0,options,LUT,L0Parabolic,R0Parabolic,POISParabolic,YOUNGParabolic,ext_verts,ext_force_status); % get initial guess using the linear solver
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
        [Ls Ltheta] = strains(rv, L0, R0);
        strainFramesl(1,:) = Ls;
        strainFramesr(1,:) = Ltheta;
    end
  
    X = fsolve(@solver_cspline_local_fast_linelastic,X0,options2,LUT,L0,R0,POIS,YOUNG,ext_verts,ext_force_status);
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)'; 
    [Ls Ltheta] = strains(rv, L0, R0);
    strainFramesl(1,:) = Ls;
    strainFramesr(1,:) = Ltheta;
    X = fsolve(@solver_cspline_local_fast_linelastic,X0,options2,LUT,L0,R0,POIS,YOUNG,ext_verts,ext_force_status);
    error = max(abs(X-X0))
    X0=X;
    frames(end+1, :) = X0;
    rv = reshape(X0',N,2)';  
    [Ls Ltheta] = strains(rv, L0, R0);
    strainFramesl(1,:) = Ls;
    strainFramesr(1,:) = Ltheta;
    rb = LUT*rv';
    D = sqrt(sum(rb.^2,2)); %Vector of edge lengths.
    rm = 0.5*(rv(2,1:end-1)+rv(2,2:end))';
    angle = atan(rb(:,2)./rb(:,1));
    
    [Tl Tr] = tensions(Ls, Ltheta, POIS, YOUNG); % column vectors
    Tl = Tl';
    Tr = Tr';
end

function [Ls Ltheta] = strains(rv, L0, R0)
    % Decide which point will separate the horizontally- and
    % vertically-oriented splines. This will happen when dr/dz = slope_cutoff.
    % slope_cutoff = -2;
    N = size(rv,2);
    slope_cutoff = -1 + 0.001;
    cut = 2;
    slope = (rv(2,3) - rv(2,1)) / (rv(1,3) - rv(1,1));
    while slope >= slope_cutoff
        cut = cut+1;
        slope = (rv(2,cut+1) - rv(2,cut-1)) / (rv(1,cut+1) - rv(1,cut-1));
    end
    [horizA horizB horizC horizD] = find_splines(rv(1,1:cut), rv(2,1:cut), 0, slope);
    [vertA vertB vertC vertD] = find_splines(-rv(2,cut:end), rv(1,cut:end), -1/slope, 0);
    arclength = zeros(N,1);
    arclength(1) = 0;
    for i = 1:cut-1
        arclength(i+1) = midpt(@(z) sqrt(1 + (3 * horizA(i) * z^2 + 2 * horizB(i) * z + horizC(i))^2), rv(1,i), rv(1,i+1), 256);
    end
    for i = cut:N-1
        arclength(i+1) = midpt(@(r) sqrt(1 + (3 * vertA(i-cut+1) * r^2 + 2 * vertB(i-cut+1) * r + vertC(i-cut+1))^2), -rv(2,i), -rv(2,i+1), 256);
    end

    S = cumsum(arclength)'; % pairs of corresponding s and s^0 coordinates
    % S0 = cumsum(L0);
    S0 = L0;
    % generate splines to interpolate the increasing function s(s^0)
    % [SA SB SC SD] = find_splines(S0, S, S(2) / S0(2), 1 ,1);
    [SA SB SC SD] = find_splines(S0, S, S(2) / S0(2), (S(end) - S(end-1)) / (S0(end) - S0(end-1)));
    [RA RB RC RD] = find_splines(S, rv(2,:), 0, 1);
    [R0A R0B R0C R0D] = find_splines(S0, R0, 0, 1);
    for i = 1:N-1
        Ls(i) = 3 * SA(i) * S0(i)^2 + 2 * SB(i) * S0(i) + SC(i); % ds/ds^0
    end
    Ls(N) = 3 * SA(end) * S0(end)^2 + 2 * SB(end) * S0(end) + SC(end);
    Ltheta(1:N-1) = rv(2,1:end-1) ./ R0(1:end-1);
    % Ltheta(N) = 1; % 0/0
    Ltheta(N) = Ls(N); % by L'Hopital's rule
end

function [Tl Tr] = tensions(Ls, Ltheta, POIS, YOUNG)
    Tl = YOUNG' ./ (1 - POIS') .* (Ls - 1);
    Tr = YOUNG' ./ (1 - POIS') .* (Ltheta - 1);
end