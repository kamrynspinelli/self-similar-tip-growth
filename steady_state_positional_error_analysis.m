function steady_state_positional_error_analysis()
    tic; % start timing execution

    %% Setup
    % parameters
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor
    initSize = 4;
    endSize = 6;
    resultsPositions = cell(endSize - initSize + 1, 4); % switch 4 to 5 when implementing degenerate parabolic case
    resultsTl = cell(endSize - initSize + 1, 4); % switch 4 to 5 when implementing degenerate parabolic case
    resultsTr = cell(endSize - initSize + 1, 4); % switch 4 to 5 when implementing degenerate parabolic case
    positionTwoNorms = zeros(endSize - initSize + 1, 4);
    positionInfNorms = zeros(endSize - initSize + 1, 4);
    TlTwoNorms = zeros(endSize - initSize + 1, 4);
    TlInfNorms = zeros(endSize - initSize + 1, 4);
    TrTwoNorms = zeros(endSize - initSize + 1, 4);
    TrInfNorms = zeros(endSize - initSize + 1, 4);
    K = 1;
    mu = 1;
    
    %% Marker point position computation
    for logSize = initSize:endSize
        [TlLinear, TrLinear, rvLinear] = run_simulation_linear(2^logSize);
        resultsPositions(logSize - initSize + 1, 1) = {rvLinear};
        resultsTl(logSize - initSize + 1, 1) = {TlLinear};
        resultsTr(logSize - initSize + 1, 1) = {TrLinear};
        [TlParabolic, TrParabolic, rvParabolic] = run_simulation_parabolic(2^logSize);
        resultsPositions(logSize - initSize + 1, 2) = {rvParabolic};
        resultsTl(logSize - initSize + 1, 2) = {TlParabolic};
        resultsTr(logSize - initSize + 1, 2) = {TrParabolic};
        [TlSpline, TrSpline, rvSpline] = run_simulation_spline(2^logSize);
        resultsPositions(logSize - initSize + 1, 3) = {rvSpline};
        resultsTl(logSize - initSize + 1, 3) = {TlSpline};
        resultsTr(logSize - initSize + 1, 3) = {TrSpline};
        [TlSplineLocal, TrSplineLocal, rvSplineLocal] = run_simulation_spline_local(2^logSize);
        resultsPositions(logSize - initSize + 1, 4) = {rvSplineLocal};
        resultsTl(logSize - initSize + 1, 4) = {TlSplineLocal};
        resultsTr(logSize - initSize + 1, 4) = {TrSplineLocal};
    end
    
    steadyRadius = 1/4*(1+sqrt(17)); % the radius of the sphere in the steady state (K = mu = P = 1)
    for logSize = initSize:endSize
        steadySphere(1,:) = steadyRadius * cos(linspace(pi/2, 0, 2^logSize +1)); % the true steady sphere on 2^logSize patches
        steadySphere(2,:) = steadyRadius * sin(linspace(pi/2, 0, 2^logSize +1));
        steadyTl = (steadyRadius^2 - 1) * ones(2^logSize, 1); % the true tensions on the sphere
        steadyTr = (steadyRadius^2 - 1) * ones(2^logSize, 1);
        positionTwoNorms(logSize - initSize + 1, 1) = norm(resultsPositions{logSize - initSize + 1, 1} - steadySphere, 2); % norms of error in the linear scheme
        positionInfNorms(logSize - initSize + 1, 1) = norm(resultsPositions{logSize - initSize + 1, 1} - steadySphere, Inf);
        TlTwoNorms(logSize - initSize + 1, 1) = norm(resultsTl{logSize - initSize + 1, 1} - steadyTl, 2);
        TlInfNorms(logSize - initSize + 1, 1) = norm(resultsTl{logSize - initSize + 1, 1} - steadyTl, Inf);
        TrTwoNorms(logSize - initSize + 1, 1) = norm(resultsTr{logSize - initSize + 1, 1} - steadyTr, 2);
        TrInfNorms(logSize - initSize + 1, 1) = norm(resultsTr{logSize - initSize + 1, 1} - steadyTr, Inf);
        positionTwoNorms(logSize - initSize + 1, 2) = norm(resultsPositions{logSize - initSize + 1, 2} - steadySphere, 2); % norms of error in the parabolic scheme
        positionInfNorms(logSize - initSize + 1, 2) = norm(resultsPositions{logSize - initSize + 1, 2} - steadySphere, Inf);
        TlTwoNorms(logSize - initSize + 1, 2) = norm(resultsTl{logSize - initSize + 1, 2} - steadyTl, 2);
        TlInfNorms(logSize - initSize + 1, 2) = norm(resultsTl{logSize - initSize + 1, 2} - steadyTl, Inf);
        TrTwoNorms(logSize - initSize + 1, 2) = norm(resultsTr{logSize - initSize + 1, 2} - steadyTr, 2);
        TrInfNorms(logSize - initSize + 1, 2) = norm(resultsTr{logSize - initSize + 1, 2} - steadyTr, Inf);
        positionTwoNorms(logSize - initSize + 1, 3) = norm(resultsPositions{logSize - initSize + 1, 3} - steadySphere, 2); % norms of error in the cubic spline scheme
        positionInfNorms(logSize - initSize + 1, 3) = norm(resultsPositions{logSize - initSize + 1, 3} - steadySphere, Inf);
        TlTwoNorms(logSize - initSize + 1, 3) = norm(resultsTl{logSize - initSize + 1, 3} - steadyTl, 2);
        TlInfNorms(logSize - initSize + 1, 3) = norm(resultsTl{logSize - initSize + 1, 3} - steadyTl, Inf);
        TrTwoNorms(logSize - initSize + 1, 3) = norm(resultsTr{logSize - initSize + 1, 3} - steadyTr, 2);
        TrInfNorms(logSize - initSize + 1, 3) = norm(resultsTr{logSize - initSize + 1, 3} - steadyTr, Inf);
        positionTwoNorms(logSize - initSize + 1, 4) = norm(resultsPositions{logSize - initSize + 1, 4} - steadySphere, 2); % norms of error in the local cubic spline scheme
        positionInfNorms(logSize - initSize + 1, 4) = norm(resultsPositions{logSize - initSize + 1, 4} - steadySphere, Inf);
%         TlTwoNorms(logSize - initSize + 1, 4) = norm(resultsTl{logSize - initSize + 1, 4} - steadyTl, 2);
%         TlInfNorms(logSize - initSize + 1, 4) = norm(resultsTl{logSize - initSize + 1, 4} - steadyTl, Inf);
%         TrTwoNorms(logSize - initSize + 1, 4) = norm(resultsTr{logSize - initSize + 1, 4} - steadyTr, 2);
%         TrInfNorms(logSize - initSize + 1, 4) = norm(resultsTr{logSize - initSize + 1, 4} - steadyTr, Inf);
        clearvars steadySphere;
    end
    
    %% Plotting
    % Plot and order of accuracy computation for marker point positions under 2-norm
    pointsArray = arrayfun(@(x) 2^x, initSize:endSize);
    fitTwoL = polyfit(log2(pointsArray), log2(positionTwoNorms(:,1)), 1);
    fitTwoP = polyfit(log2(pointsArray), log2(positionTwoNorms(:,2)), 1);
    fitTwoS = polyfit(log2(pointsArray), log2(positionTwoNorms(:,3)), 1);
    fitTwoSL = polyfit(log2(pointsArray), log2(positionTwoNorms(:,4)), 1);
    orderTwoL = -fitTwoL(1);
    orderTwoP = -fitTwoP(1);
    orderTwoS = -fitTwoS(1);
    orderTwoSL = -fitTwoSL(1);
    hold on;
    scatter(pointsArray, positionTwoNorms(:,1), 'filled');
    scatter(pointsArray, positionTwoNorms(:,2), 'filled');
    scatter(pointsArray, positionTwoNorms(:,3), 'filled');
    scatter(pointsArray, positionTwoNorms(:,4), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoL(2) * x^fitTwoL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoP(2) * x^fitTwoP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoS(2) * x^fitTwoS(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoSL(2) * x^fitTwoSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', 'Cubic splines', 'Cubic splines (local)', ...
        ['Linear fitting: slope = ', num2str(fitTwoL(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoP(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoS(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoSL(1))]);
    title('Computed positional error (2-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_position_error_analysis_2norm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for marker point positions in the linear model under 2-norm is: %f.\n", orderTwoL);
    fprintf("The computed order of accuracy for marker point positions in the parabolic model under 2-norm is: %f.\n", orderTwoP);
    fprintf("The computed order of accuracy for marker point positions in the cubic spline model under 2-norm is: %f.\n", orderTwoS);
    fprintf("The computed order of accuracy for marker point positions in the local cubic spline model under 2-norm is: %f.\n", orderTwoSL);
    
    % Plot and order of accuracy computation for marker point positions under infinity-norm
    fitInfL = polyfit(log2(pointsArray), log2(positionInfNorms(:,1)), 1);
    fitInfP = polyfit(log2(pointsArray), log2(positionInfNorms(:,2)), 1);
    fitInfS = polyfit(log2(pointsArray), log2(positionInfNorms(:,3)), 1);
    orderInfL = -fitInfL(1);
    orderInfP = -fitInfP(1);
    orderInfS = -fitInfS(1);
    hold on;
    scatter(pointsArray, positionInfNorms(:,1), 'filled');
    scatter(pointsArray, positionInfNorms(:,2), 'filled');
    scatter(pointsArray, positionInfNorms(:,3), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfL(2) * x^fitInfL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfP(2) * x^fitInfP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfS(2) * x^fitInfS(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', 'Cubic splines', ...
        ['Linear fitting: slope = ', num2str(fitInfL(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfP(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfS(1))]);
    title('Computed positional error (infinity-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_position_error_analysis_Infnorm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for marker point positions in the linear model under infinity-norm is: %f.\n", orderInfL);
    fprintf("The computed order of accuracy for marker point positions in the parabolic model under infinity-norm is: %f.\n", orderInfP);
    fprintf("The computed order of accuracy for marker point positions in the cubic spline model under infinity-norm is: %f.\n", orderInfS);
    
    % Plot and order of accuracy computation for meridional tension under 2-norm
    pointsArray = arrayfun(@(x) 2^x, initSize:endSize);
    fitTwoL = polyfit(log2(pointsArray), log2(TlTwoNorms(:,1)), 1);
    fitTwoP = polyfit(log2(pointsArray), log2(TlTwoNorms(:,2)), 1);
    orderTwoL = -fitTwoL(1);
    orderTwoP = -fitTwoP(1);
    hold on;
    scatter(pointsArray, TlTwoNorms(:,1), 'filled');
    scatter(pointsArray, TlTwoNorms(:,2), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoL(2) * x^fitTwoL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoP(2) * x^fitTwoP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', ...
        ['Linear fitting: slope = ', num2str(fitTwoL(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoP(1))]);
    title('Computed meridional tension error (2-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_Tl_error_analysis_2norm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for meridional tensions in the linear scheme under 2-norm is: %f.\n", orderTwoL);
    fprintf("The computed order of accuracy for meridional tensions in the parabolic scheme under 2-norm is: %f.\n", orderTwoP);
    
    % Plot and order of accuracy computation for meridional tension under infinity-norm
    fitInfL = polyfit(log2(pointsArray), log2(TlInfNorms(:,1)), 1);
    fitInfP = polyfit(log2(pointsArray), log2(TlInfNorms(:,2)), 1);
    orderInfL = -fitInfL(1);
    orderInfP = -fitInfP(1);
    hold on;
    scatter(pointsArray, TlInfNorms(:,1), 'filled');
    scatter(pointsArray, TlInfNorms(:,2), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfL(2) * x^fitInfL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfP(2) * x^fitInfP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', ...
        ['Linear fitting: slope = ', num2str(fitInfL(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfP(1))]);
    title('Computed meridional tension error (infinity-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_Tl_error_analysis_Infnorm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for meridional tensions in the linear scheme under infinity-norm is: %f.\n", orderInfL);
    fprintf("The computed order of accuracy for meridional tensions in the parabolic scheme under infinity-norm is: %f.\n", orderInfP);
    
    % Plot and order of accuracy computation for circumferential tension under 2-norm
    pointsArray = arrayfun(@(x) 2^x, initSize:endSize);
    fitTwoL = polyfit(log2(pointsArray), log2(TrTwoNorms(:,1)), 1);
    fitTwoP = polyfit(log2(pointsArray), log2(TrTwoNorms(:,2)), 1);
    orderTwoL = -fitTwoL(1);
    orderTwoP = -fitTwoP(1);
    hold on;
    scatter(pointsArray, TrTwoNorms(:,1), 'filled');
    scatter(pointsArray, TrTwoNorms(:,2), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoL(2) * x^fitTwoL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoP(2) * x^fitTwoP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', ...
        ['Linear fitting: slope = ', num2str(fitTwoL(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoP(1))]);
    title('Computed circumferential tension error (2-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_Tr_error_analysis_2norm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for circumferential tensions in the linear scheme under 2-norm is: %f.\n", orderTwoL);
    fprintf("The computed order of accuracy for circumferential tensions in the parabolic scheme under 2-norm is: %f.\n", orderTwoP);
    
    % Plot and order of accuracy computation for circumferential tension under infinity-norm
    fitInfL = polyfit(log2(pointsArray), log2(TrInfNorms(:,1)), 1);
    fitInfP = polyfit(log2(pointsArray), log2(TrInfNorms(:,2)), 1);
    orderInfL = -fitInfL(1);
    orderInfP = -fitInfP(1);
    hold on;
    scatter(pointsArray, TrInfNorms(:,1), 'filled');
    scatter(pointsArray, TrInfNorms(:,2), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfL(2) * x^fitInfL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfP(2) * x^fitInfP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', ...
        ['Linear fitting: slope = ', num2str(fitInfL(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfP(1))]);
    title('Computed circumferential tension error (infinity-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_Tr_error_analysis_Infnorm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of accuracy for circumferential tensions in the linear scheme under infinity-norm is: %f.\n", orderInfL);
    fprintf("The computed order of accuracy for circumferential tensions in the parabolic scheme under infinity-norm is: %f.\n", orderInfP);
    
    toc % show the execution time
end

function [Tl, Tr, rv] = run_simulation_parabolic(N)
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor
    
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
    Inc =  1;

    % external force parameters
    external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
    Fext = 0;
    variance_factor = 1;

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
    
    rv = 2 * rv; % for the initial guess to ode45, pick an expanded shell
    % run the simulation
    [ Tl,Tr,rv,error ] = equilibrium_parabolic_imperative(rv,LUT,L0Parabolic,R0Parabolic,K,MU,ext_verts,0);
end

function [Tl, Tr, rv] = run_simulation_linear(N)
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor
    
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
    Inc =  1;

    % external force parameters
    external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
    Fext = 0;
    variance_factor = 1;

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
    end
    ext_verts = find(sum(adj,1)==1);
    
    rv = 2 * rv; % for the initial guess to ode45, pick an expanded shell
    % run the simulation
    [ Tl,Tr,rv,error ] = equilibrium_linear(rv,LUT,L0Linear,R0Linear,K,MU,ext_verts,0); % linear segments
end

function [Tl, Tr, rv] = run_simulation_spline(N)
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor
    
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
    Inc =  1;

    % external force parameters
    external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
    Fext = 0;
    variance_factor = 1;

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
    end
    [L0Spline R0Spline] = spline_intrinsics(rv);
    ext_verts = find(sum(adj,1)==1);
    
    rv = 2 * rv; % for the initial guess to ode45, pick an expanded shell
    % run the simulation
    [ Tl,Tr,rv,error ] = equilibrium_cspline(rv,LUT,L0Spline,R0Spline,K,MU,ext_verts,0); % linear segments
end

function [Tl, Tr, rv] = run_simulation_spline_local(N)
    global Tol Rtol TolFun TolX Inc L P
    global Fext external_force_mode variance_factor
    
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
    Inc = 0.05;

    % external force parameters
    external_force_mode = "gaussian"; % one of "tip_only" or "gaussian"
    Fext = 0;
    variance_factor = 1;

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
    [ Tl,Tr,rv,error ] = equilibrium_cspline_local(rv,LUT,L0SplineLocal,R0SplineLocal,KSplineLocal,MUSplineLocal,ext_verts,0,L0Parabolic,R0Parabolic,K,MU);
end