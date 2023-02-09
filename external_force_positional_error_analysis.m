function external_force_positional_error_analysis(force_mode, force_magnitude, force_dist)
    % force_mode: one of "tip_only", "gaussian", "triangular"
    % Fext: the magnitude of force to be applied
    % dist: the distributional parameter of the force application (variance
    % for Gaussian distribution, etc.)

    tic; % start timing execution

    %% Setup
    % parameters
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor triangle_half_width
    initSize = 4;
    endSize = 6;
    resultsPositions = cell(endSize - initSize + 1, 4);
    resultsTl = cell(endSize - initSize + 1, 4);
    resultsTr = cell(endSize - initSize + 1, 4);
    positionTwoNorms = zeros(endSize - initSize + 1, 4);
    positionInfNorms = zeros(endSize - initSize + 1, 4);
    TlTwoNorms = zeros(endSize - initSize + 1, 4);
    TlInfNorms = zeros(endSize - initSize + 1, 4);
    TrTwoNorms = zeros(endSize - initSize + 1, 4);
    TrInfNorms = zeros(endSize - initSize + 1, 4);
    K = 1;
    mu = 1;
    
    external_force_mode = force_mode;
    Fext = force_magnitude;
    if force_mode == "gaussian"
        variance_factor = force_dist;
    elseif force_mode == "triangular"
        triangle_half_width = force_dist;
    end
    
    %% Marker point position computation
    for logSize = initSize:endSize
%         [TlLinear, TrLinear, rvLinear] = run_simulation_linear(2^logSize);
%         resultsPositions(logSize - initSize + 1, 1) = {rvLinear};
%         resultsTl(logSize - initSize + 1, 1) = {TlLinear};
%         resultsTr(logSize - initSize + 1, 1) = {TrLinear};
        [TlParabolic, TrParabolic, rvParabolic] = run_simulation_parabolic(2^logSize);
        resultsPositions(logSize - initSize + 1, 2) = {rvParabolic};
        resultsTl(logSize - initSize + 1, 2) = {TlParabolic};
        resultsTr(logSize - initSize + 1, 2) = {TrParabolic};
%         [TlParabolicDegenerate, TrParabolicDegenerate, rvParabolicDegenerate] = run_simulation_parabolic_degenerate(2^logSize, a, b);
%         resultsPositions(logSize - initSize + 1, 3) = {rvParabolicDegenerate};
%         resultsTl(logSize - initSize + 1, 3) = {TlParabolicDegenerate};
%         resultsTr(logSize - initSize + 1, 3) = {TrParabolicDegenerate};
        [TlSplineLocal, TrSplineLocal, rvSplineLocal] = run_simulation_spline_local(2^logSize);
        resultsPositions(logSize - initSize + 1, 4) = {rvSplineLocal};
        resultsTl(logSize - initSize + 1, 4) = {TlSplineLocal};
        resultsTr(logSize - initSize + 1, 4) = {TrSplineLocal};
    end
    
    for logSize = initSize:endSize
%         positionTwoNorms(logSize - initSize + 1, 1) = norm(resultsPositions{logSize - initSize + 1, 1} - ...
%             trim_to(resultsPositions{endSize - initSize + 1, 1}, endSize, logSize), 2); % norms of error in the linear scheme
%         positionInfNorms(logSize - initSize + 1, 1) = norm(resultsPositions{logSize - initSize + 1, 1} - ...
%             trim_to(resultsPositions{endSize - initSize + 1, 1}, endSize, logSize), Inf);
%         TlTwoNorms(logSize - initSize + 1, 1) = norm(resultsTl{logSize - initSize + 1, 1} - ...
%             average_to(resultsTl{endSize - initSize + 1, 1}, endSize, logSize), 2);
%         TlInfNorms(logSize - initSize + 1, 1) = norm(resultsTl{logSize - initSize + 1, 1} - ...
%             average_to(resultsTl{endSize - initSize + 1, 1}, endSize, logSize), Inf);
%         TrTwoNorms(logSize - initSize + 1, 1) = norm(resultsTr{logSize - initSize + 1, 1} - ...
%             average_to(resultsTr{endSize - initSize + 1, 1}, endSize, logSize), 2);
%         TrInfNorms(logSize - initSize + 1, 1) = norm(resultsTr{logSize - initSize + 1, 1} - ...
%             average_to(resultsTr{endSize - initSize + 1, 1}, endSize, logSize), Inf);
        positionTwoNorms(logSize - initSize + 1, 2) = norm(resultsPositions{logSize - initSize + 1, 2} - ...
            trim_to(resultsPositions{endSize - initSize + 1, 2}, endSize, logSize), 2); % norms of error in the parabolic scheme
        positionInfNorms(logSize - initSize + 1, 2) = norm(resultsPositions{logSize - initSize + 1, 2} - ...
            trim_to(resultsPositions{endSize - initSize + 1, 2}, endSize, logSize), Inf);
        TlTwoNorms(logSize - initSize + 1, 2) = norm(resultsTl{logSize - initSize + 1, 2} - ...
            average_to(resultsTl{endSize - initSize + 1, 2}, endSize, logSize), 2);
        TlInfNorms(logSize - initSize + 1, 2) = norm(resultsTl{logSize - initSize + 1, 2} - ...
            average_to(resultsTl{endSize - initSize + 1, 2}, endSize, logSize), Inf);
        TrTwoNorms(logSize - initSize + 1, 2) = norm(resultsTr{logSize - initSize + 1, 2} - ...
            average_to(resultsTr{endSize - initSize + 1, 2}, endSize, logSize), 2);
        TrInfNorms(logSize - initSize + 1, 2) = norm(resultsTr{logSize - initSize + 1, 2} - ...
            average_to(resultsTr{endSize - initSize + 1, 2}, endSize, logSize), Inf);
%         positionTwoNorms(logSize - initSize + 1, 3) = norm(resultsPositions{logSize - initSize + 1, 3} - ...
%             trim_to(resultsPositions{endSize - initSize + 1, 3}, endSize, logSize), 2); % norms of error in the degenerate parabolic scheme
%         positionInfNorms(logSize - initSize + 1, 3) = norm(resultsPositions{logSize - initSize + 1, 3} - ...
%             trim_to(resultsPositions{endSize - initSize + 1, 3}, endSize, logSize), Inf);
%         TlTwoNorms(logSize - initSize + 1, 3) = norm(resultsTl{logSize - initSize + 1, 3} - ...
%             average_to(resultsTl{endSize - initSize + 1, 3}, endSize, logSize), 2);
%         TlInfNorms(logSize - initSize + 1, 3) = norm(resultsTl{logSize - initSize + 1, 3} - ...
%             average_to(resultsTl{endSize - initSize + 1, 3}, endSize, logSize), Inf);
%         TrTwoNorms(logSize - initSize + 1, 3) = norm(resultsTr{logSize - initSize + 1, 3} - ...
%             average_to(resultsTr{endSize - initSize + 1, 3}, endSize, logSize), 2);
%         TrInfNorms(logSize - initSize + 1, 3) = norm(resultsTr{logSize - initSize + 1, 3} - ...
%             average_to(resultsTr{endSize - initSize + 1, 3}, endSize, logSize), Inf);
        positionTwoNorms(logSize - initSize + 1, 4) = norm(resultsPositions{logSize - initSize + 1, 4} - ...
            trim_to(resultsPositions{endSize - initSize + 1, 4}, endSize, logSize), 2); % norms of error in the local cubic spline scheme
        positionInfNorms(logSize - initSize + 1, 4) = norm(resultsPositions{logSize - initSize + 1, 4} - ...
            trim_to(resultsPositions{endSize - initSize + 1, 4}, endSize, logSize), Inf);
%         TlTwoNorms(logSize - initSize + 1, 4) = norm(resultsTl{logSize - initSize + 1, 4} - ...
%             average_to(resultsTl{endSize - initSize + 1, 4}, endSize, logSize), 2);
%         TlInfNorms(logSize - initSize + 1, 4) = norm(resultsTl{logSize - initSize + 1, 4} - ...
%             average_to(resultsTl{endSize - initSize + 1, 4}, endSize, logSize), Inf);
%         TrTwoNorms(logSize - initSize + 1, 4) = norm(resultsTr{logSize - initSize + 1, 4} - ...
%             average_to(resultsTr{endSize - initSize + 1, 4}, endSize, logSize), 2);
%         TrInfNorms(logSize - initSize + 1, 4) = norm(resultsTr{logSize - initSize + 1, 4} - ...
%             average_to(resultsTr{endSize - initSize + 1, 4}, endSize, logSize), Inf);
    end
    
    %% Plotting
    % Plot and order of accuracy computation for marker point positions under 2-norm
    pointsArray = arrayfun(@(x) 2^x, initSize:endSize-1);
%     fitTwoL = polyfit(log2(pointsArray), log2(positionTwoNorms(1:end-1,1)), 1);
    fitTwoP = polyfit(log2(pointsArray), log2(positionTwoNorms(1:end-1,2)), 1);
%     fitTwoPD = polyfit(log2(pointsArray), log2(positionTwoNorms(1:end-1,3)), 1);
    fitTwoSL = polyfit(log2(pointsArray), log2(positionTwoNorms(1:end-1,4)), 1);
%     orderTwoL = -fitTwoL(1);
    orderTwoP = -fitTwoP(1);
%     orderTwoPD = -fitTwoPD(1);
    orderTwoSL = -fitTwoSL(1);
    hold on;
%     scatter(pointsArray, positionTwoNorms(1:end-1,1), 'filled');
    scatter(pointsArray, positionTwoNorms(1:end-1,2), 'filled');
%     scatter(pointsArray, positionTwoNorms(1:end-1,3), 'filled');
    scatter(pointsArray, positionTwoNorms(1:end-1,4), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoP(2) * x^fitTwoP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoSL(2) * x^fitTwoSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTwoL(2) * x^fitTwoL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTwoPD(2) * x^fitTwoPD(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Parabolic arcs', 'Cubic splines (local)', ...
        ['Linear fitting: slope = ', num2str(fitTwoP(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoSL(1))]);
%         ['Linear fitting: slope = ', num2str(fitTwoL(1))], ...
%         ['Linear fitting: slope = ', num2str(fitTwoPD(1))], ...
    if external_force_mode == "gaussian"
        % title(['Computed positional error (2-norm) with external force (Fext=', num2str(force_magnitude), ',var=', num2str(variance_factor), ')']);
    elseif external_force_mode == "triangular"
        % title(['Computed positional error (2-norm) with external force (Fext=', num2str(force_magnitude), ',width=', num2str(triangle_half_width), ')']);
    end
    pbaspect([1 1 1]);
    if external_force_mode == "gaussian"
        exportgraphics(gcf, ['media/external_force_mode_gaussian_', 'Fext_', num2str(force_magnitude), '_var_', num2str(variance_factor), '_position_error_analysis_2norm.png'], 'Resolution', 600);
    elseif external_force_mode == "triangular"
        exportgraphics(gcf, ['media/external_force_mode_triangular_', 'Fext_', num2str(force_magnitude), '_width_', num2str(triangle_half_width), '_position_error_analysis_2norm.png'], 'Resolution', 600);
    end
    close all;
    
%     fprintf("The computed order of accuracy for marker point positions in the linear model under 2-norm is: %f.\n", orderTwoL);
    fprintf("The computed order of accuracy for marker point positions in the parabolic model under 2-norm is: %f.\n", orderTwoP);
%     fprintf("The computed order of accuracy for marker point positions in the degenerate parabolic model under 2-norm is: %f.\n", orderTwoPD);
    fprintf("The computed order of accuracy for marker point positions in the (local) cubic spline model under 2-norm is: %f.\n", orderTwoSL);
    
    % Plot and order of accuracy computation for marker point positions under infinity-norm
%     fitInfL = polyfit(log2(pointsArray), log2(positionInfNorms(1:end-1,1)), 1);
    fitInfP = polyfit(log2(pointsArray), log2(positionInfNorms(1:end-1,2)), 1);
    % fitInfPD = polyfit(log2(pointsArray), log2(positionInfNorms(1:end-1,3)), 1);
    fitInfSL = polyfit(log2(pointsArray), log2(positionInfNorms(1:end-1,4)), 1);
%     orderInfL = -fitInfL(1);
    orderInfP = -fitInfP(1);
    % orderInfPD = -fitInfPD(1);
    orderInfSL = -fitInfSL(1);
    hold on;
%     scatter(pointsArray, positionInfNorms(1:end-1,1), 'filled');
    scatter(pointsArray, positionInfNorms(1:end-1,2), 'filled');
    % scatter(pointsArray, positionInfNorms(1:end-1,3), 'filled');
    scatter(pointsArray, positionInfNorms(1:end-1,4), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfP(2) * x^fitInfP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfSL(2) * x^fitInfSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitInfL(2) * x^fitInfL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...        
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitInfPD(2) * x^fitInfPD(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Parabolic arcs', 'Cubic splines (local)', ...
        ['Linear fitting: slope = ', num2str(fitInfP(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfSL(1))]);
%         ['Linear fitting: slope = ', num2str(fitInfL(1))], ...
    %     ['Linear fitting: slope = ', num2str(fitInfPD(1))]
    if external_force_mode == "gaussian"
        % title(['Computed positional error (\infty-norm) with external force (Fext=', num2str(force_magnitude), ',var=', num2str(variance_factor), ')']);
    elseif external_force_mode == "triangular"
        % title(['Computed positional error (\infty norm) with external force (Fext=', num2str(force_magnitude), ',width=', num2str(triangle_half_width), ')']);
    end
    pbaspect([1 1 1]);
    if external_force_mode == "gaussian"
        exportgraphics(gcf, ['media/external_force_mode_gaussian_', 'Fext_', num2str(force_magnitude), '_var_', num2str(variance_factor), '_position_error_analysis_Infnorm.png'], 'Resolution', 600);
    elseif external_force_mode == "triangular"
        exportgraphics(gcf, ['media/external_force_mode_triangular_', 'Fext_', num2str(force_magnitude), '_width_', num2str(triangle_half_width), '_position_error_analysis_Infnorm.png'], 'Resolution', 600);
    end
    close all;
    
%     fprintf("The computed order of accuracy for marker point positions in the linear model under infinity-norm is: %f.\n", orderInfL);
    fprintf("The computed order of accuracy for marker point positions in the parabolic model under infinity-norm is: %f.\n", orderInfP);
    % fprintf("The computed order of accuracy for marker point positions in the degenerate parabolic model under infinity-norm is: %f.\n", orderInfPD);
    fprintf("The computed order of accuracy for marker point positions in the (local) cubic splines model under infinity-norm is: %f.\n", orderInfSL);
end

function [Tl, Tr, rv] = run_simulation_parabolic(N)
    global Tol Rtol TolFun TolX Inc L P
    
    % dimensional parameters
    L = 1;%characteristic length
    P = 1;%characteristic pressure

    % Bulk Modulus
    k = 1;
    % Shear Modulus
    mu = 1;
    % Uniform material property
    K = k*ones(N,1);
    MU = mu*ones(N,1);

    % computational parameters
    Rtol = 1e-13;
    Tol =  1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    % Inc =  1;
    Inc = 0.25;

    % initialize cells
    R = 1;%the intrinsic radius, if a sphere is to be generated
    L0Linear = ones(N,1);
    R0Linear = ones(N,1);
    L0Parabolic = ones(N,1);
    R0Parabolic = ones(N,1);
    [ rv, adj ] = generate_sphere(R,N);
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
       L0Parabolic(i) = initialArcs(i+1).arclength; % parabolic arcs
       R0Parabolic(i) = initialArcs(i+1).vert(2);
    end
    ext_verts = find(sum(adj,1)==1);
    
    rv = 2 * rv;
    % run the simulation
    [ Tl,Tr,rv,error ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0);
    % now add force
    [ Tl,Tr,rv,error ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,1);
end

function [Tl, Tr, rv] = run_simulation_linear(N)
    global Tol Rtol TolFun TolX Inc L P
    
    % dimensional parameters
    L = 1;%characteristic length
    P = 1;%characteristic pressure

    % Bulk Modulus
    k = 1;
    % Shear Modulus
    mu = 1;
    % Uniform material property
    K = k*ones(N,1);
    MU = mu*ones(N,1);

    % computational parameters
    Rtol = 1e-13;
    Tol =  1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    % Inc =  1;
    Inc = 0.25;

    % initialize cells
    R = 1;%the intrinsic radius, if a sphere is to be generated
    L0Linear = ones(N,1);
    R0Linear = ones(N,1);
    L0Parabolic = ones(N,1);
    R0Parabolic = ones(N,1);
    [ rv, adj ] = generate_sphere(R,N);
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
    end
    ext_verts = find(sum(adj,1)==1);
    
    rv = 2 * rv;
    % run the simulation
    [ Tl,Tr,rv,error ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
    % now add force
    [ Tl,Tr,rv,error ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,1); % linear segments
end

function [Tl, Tr, rv] = run_simulation_spline_local(N)
    global Tol Rtol TolFun TolX Inc L P
    
    % dimensional parameters
    L = 1;%characteristic length
    P = 1;%characteristic pressure

    % Bulk Modulus
    k = 1;
    % Shear Modulus
    mu = 1;
    % Uniform material property
    K = k*ones(N,1);
    MU = mu*ones(N,1);

    % computational parameters
    Rtol = 1e-13;
    Tol =  1e-3;
    TolFun = 1e-16;
    TolX = 1e-16;
    % Inc = 1;
    Inc =  0.25;
    % Inc = 0.1;
    % Inc = 0.05;
    % Inc = 0.03;

    % initialize cells
    R = 1;%the intrinsic radius, if a sphere is to be generated
    L0Linear = ones(N,1);
    R0Linear = ones(N,1);
    L0Parabolic = ones(N,1);
    R0Parabolic = ones(N,1);
    [ rv, adj ] = generate_sphere(R,N);
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
    ext_verts = find(sum(adj,1)==1);
    
    % rv = 2 * rv; % for the initial guess to ode45, pick an expanded shell
    % run the simulation
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU);
    % now add force
    [ Tl,Tr,rv,error,frames1,strainFramesl1,strainFramesr1 ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,1,L0Parabolic,R0Parabolic,K,MU);
end

function trimmed = trim_to(mat, fromSize, toSize)
    % trims a k-by-2^(fromSize)+1 matrix to a k-by-2^(toSize)+1 matrix by
    % retaining only columns of the form 1 + i * 2^(fromSize - toSize)
    trimmed = mat(:, 1:2^(fromSize - toSize):(2^fromSize+1));
end

function avgs = average_to(vec, fromSize, toSize)
    % averages a 2^(fromSize)-vector to a 2^(toSize)-vector by computing 
    % the average of the elements from 1 to 2^(fromSize-toSize), the 
    % average of the elements from 2^(fromSize-toSize) to 2 * 2^(fromSize-toSize),
    % and so forth. we assume the vectors to be column vectors.
    avgs = zeros(2^toSize, 1);
    for i = 1:2^toSize
        avgs(i) = mean(vec(((i-1) * 2^(fromSize - toSize) + 1):(i * 2^(fromSize - toSize))));
    end
end