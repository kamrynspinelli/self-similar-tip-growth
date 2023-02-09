function steady_state_force_error_analysis()
    %% Setup
    initSize = 4;
    endSize = 10;
    results = cell(endSize - initSize + 1, 5);
    overallTwoNorms = zeros(endSize - initSize + 1, 5);
    overallInfNorms = zeros(endSize - initSize + 1, 5);
    tipNorms = zeros(endSize - initSize + 1, 5);
    K = 1;
    mu = 1;
    
    %% Force computation
    for logSize = initSize:endSize
        [fL, fP, fPD, fS, fSL] = steady_state_forces(2^logSize, K, mu);
        results(logSize - initSize + 1, 1) = {fL};
        results(logSize - initSize + 1, 2) = {fP};
        results(logSize - initSize + 1, 3) = {fPD};
        results(logSize - initSize + 1, 4) = {fS};
        results(logSize - initSize + 1, 5) = {fSL};
    end
    
    for logSize = initSize:endSize
        fL = reshape(results{logSize - initSize + 1, 1}, 2 * 2^logSize + 2, 1);
        overallTwoNorms(logSize - initSize + 1, 1) = norm(fL(2:end-1), 2);
        overallInfNorms(logSize - initSize + 1, 1) = norm(fL(2:end-1), Inf);
        tipNorms(logSize - initSize + 1, 1) = norm(fL(end-1)); % the z-component of the tip force
        fP = reshape(results{logSize - initSize + 1, 2}, 2 * 2^logSize + 2, 1);
        overallTwoNorms(logSize - initSize + 1, 2) = norm(fP(2:end-1), 2);
        overallInfNorms(logSize - initSize + 1, 2) = norm(fP(2:end-1), Inf);
        tipNorms(logSize - initSize + 1, 2) = norm(fP(end-1)); % the z-component of the tip force
        fPD = reshape(results{logSize - initSize + 1, 3}, 2 * 2^logSize + 2, 1);
        overallTwoNorms(logSize - initSize + 1, 3) = norm(fPD(2:end-1), 2);
        overallInfNorms(logSize - initSize + 1, 3) = norm(fPD(2:end-1), Inf);
        tipNorms(logSize - initSize + 1, 3) = norm(fPD(end-1)); % the z-component of the tip force
        fS = reshape(results{logSize - initSize + 1, 4}, 2 * 2^logSize + 2, 1);
        overallTwoNorms(logSize - initSize + 1, 4) = norm(fS(2:end-1), 2);
        overallInfNorms(logSize - initSize + 1, 4) = norm(fS(2:end-1), Inf);
        tipNorms(logSize - initSize + 1, 4) = norm(fS(end-1)); % the z-component of the tip force
        fSL = reshape(results{logSize - initSize + 1, 5}, 2 * 2^logSize + 2, 1);
        overallTwoNorms(logSize - initSize + 1, 5) = norm(fSL(2:end-1), 2);
        overallInfNorms(logSize - initSize + 1, 5) = norm(fSL(2:end-1), Inf);
        tipNorms(logSize - initSize + 1, 5) = norm(fSL(end-1)); % the z-component of the tip force
    end
    
    %% Plotting
    % Plot and order of accuracy computation for overall force under 2-norm
    pointsArray = arrayfun(@(x) 2^x, initSize:endSize);
    fitTwoL = polyfit(log2(pointsArray), log2(overallTwoNorms(:,1)), 1);
    fitTwoP = polyfit(log2(pointsArray), log2(overallTwoNorms(:,2)), 1);
%     fitTwoPD = polyfit(log2(pointsArray), log2(overallTwoNorms(:,3)), 1);
%     fitTwoS = polyfit(log2(pointsArray), log2(overallTwoNorms(:,4)), 1);
    fitTwoSL = polyfit(log2(pointsArray), log2(overallTwoNorms(:,5)), 1);
    orderTwoL = -fitTwoL(1);
    orderTwoP = -fitTwoP(1);
%     orderTwoPD = -fitTwoPD(1);
%     orderTwoS = -fitTwoS(1);
    orderTwoSL = -fitTwoSL(1);
    hold on;
    scatter(pointsArray, overallTwoNorms(:,1), 'filled');
    scatter(pointsArray, overallTwoNorms(:,2), 'filled');
%     scatter(pointsArray, overallTwoNorms(:,3), 'filled');
%     scatter(pointsArray, overallTwoNorms(:,4), 'filled');
    scatter(pointsArray, overallTwoNorms(:,5), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoL(2) * x^fitTwoL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoP(2) * x^fitTwoP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTwoSL(2) * x^fitTwoSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTwoPD(2) * x^fitTwoPD(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTwoS(2) * x^fitTwoS(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', 'Cubic splines (local)', ...
        ['Linear fitting: slope = ', num2str(fitTwoL(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoP(1))], ...
        ['Linear fitting: slope = ', num2str(fitTwoSL(1))]);
%         ['Linear fitting: slope = ', num2str(fitTwoPD(1))], ...
%         ['Linear fitting: slope = ', num2str(fitTwoS(1))], ...
    title('Computed overall force (2-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_overall_force_error_analysis_2norm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of convergence for overall forces in the linear scheme under 2-norm is: %f.\n", orderTwoL);
    fprintf("The computed order of convergence for overall forces in the parabolic scheme under 2-norm is: %f.\n", orderTwoP);
%     fprintf("The computed order of convergence for overall forces in the degenerate parabolic scheme under 2-norm is: %f.\n", orderTwoPD);
%     fprintf("The computed order of convergence for overall forces in the cubic spline scheme under 2-norm is: %f.\n", orderTwoS);
    fprintf("The computed order of convergence for overall forces in the cubic spline (local) scheme under 2-norm is: %f.\n", orderTwoSL);
    
    % Plot and order of accuracy computation for overall force under infinity-norm
    fitInfL = polyfit(log2(pointsArray), log2(overallInfNorms(:,1)), 1);
    fitInfP = polyfit(log2(pointsArray), log2(overallInfNorms(:,2)), 1);
%     fitInfPD = polyfit(log2(pointsArray), log2(overallInfNorms(:,3)), 1);
    fitInfSL = polyfit(log2(pointsArray), log2(overallInfNorms(:,5)), 1);
    orderInfL = -fitInfL(1);
    orderInfP = -fitInfP(1);
%     orderInfPD = -fitInfPD(1);
    orderInfSL = -fitInfSL(1);
    hold on;
    scatter(pointsArray, overallInfNorms(:,1), 'filled');
    scatter(pointsArray, overallInfNorms(:,2), 'filled');
%     scatter(pointsArray, overallInfNorms(:,3), 'filled');
    scatter(pointsArray, overallInfNorms(:,5), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfL(2) * x^fitInfL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfP(2) * x^fitInfP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitInfSL(2) * x^fitInfSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', 'Cubic splines', ...
        ['Linear fitting: slope = ', num2str(fitInfL(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfP(1))], ...
        ['Linear fitting: slope = ', num2str(fitInfSL(1))]);
    title('Computed overall force (infinity-norm) on the sphere');
    exportgraphics(gcf, 'media/sphere_overall_force_error_analysis_Infnorm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of convergence for overall forces in the linear scheme under infinity-norm is: %f.\n", orderInfL);
    fprintf("The computed order of convergence for overall forces in the parabolic scheme under infinity-norm is: %f.\n", orderInfP);
%     fprintf("The computed order of convergence for overall forces in the degenerate parabolic scheme under infinity-norm is: %f.\n", orderInfPD);
    fprintf("The computed order of convergence for overall forces in the cubic spline (local) scheme under infinity-norm is: %f.\n", orderInfSL);
    
    % Plot and order of accuracy computation for tip force
    fitTipL = polyfit(log2(pointsArray), log2(tipNorms(:,1)), 1);
    fitTipP = polyfit(log2(pointsArray), log2(tipNorms(:,2)), 1);
%     fitTipPD = polyfit(log2(pointsArray), log2(tipNorms(:,3)), 1);
%     fitTipS = polyfit(log2(pointsArray), log2(tipNorms(:,4)), 1);
    fitTipSL = polyfit(log2(pointsArray), log2(tipNorms(:,5)), 1);
    orderTipL = -fitTipL(1);
    orderTipP = -fitTipP(1);
%     orderTipPD = -fitTipPD(1);
%     orderTipS = -fitTipS(1);
    orderTipSL = -fitTipSL(1);
    hold on;
    scatter(pointsArray, tipNorms(:,1), 'filled');
    scatter(pointsArray, tipNorms(:,2), 'filled');
%     scatter(pointsArray, tipNorms(:,3), 'filled');
%     scatter(pointsArray, tipNorms(:,4), 'filled');
    scatter(pointsArray, tipNorms(:,5), 'filled');
    plot(linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTipL(2) * x^fitTipL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTipP(2) * x^fitTipP(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        linspace(pointsArray(1), pointsArray(end), 64), ...
        arrayfun(@(x) 2^fitTipSL(2) * x^fitTipSL(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
        'LineWidth', 2.0);
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTipPD(2) * x^fitTipPD(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
%         linspace(pointsArray(1), pointsArray(end), 64), ...
%         arrayfun(@(x) 2^fitTipS(2) * x^fitTipS(1), linspace(pointsArray(1), pointsArray(end), 64)), ...
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    axis([pointsArray(1) pointsArray(end) 0 1]); axis('auto y');
    legend('Linear segments', 'Parabolic arcs', 'Cubic splines (local)', ...
        ['Linear fitting: slope = ', num2str(fitTipL(1))], ...
        ['Linear fitting: slope = ', num2str(fitTipP(1))], ...
        ['Linear fitting: slope = ', num2str(fitTipSL(1))]);
%         ['Linear fitting: slope = ', num2str(fitTipPD(1))], ...
%         ['Linear fitting: slope = ', num2str(fitTipS(1))], ...
    title('Computed tip force on the sphere');
    exportgraphics(gcf, 'media/sphere_tip_force_error_analysis_Infnorm.png', 'Resolution', 600);
    close all;
    
    fprintf("The computed order of convergence for tip force in the linear scheme is: %f.\n", orderTipL);
    fprintf("The computed order of convergence for tip force in the parabolic scheme is: %f.\n", orderTipP);
%     fprintf("The computed order of convergence for tip force in the degenerate parabolic scheme is: %f.\n", orderTipPD);
%     fprintf("The computed order of convergence for tip force in the cubic spline scheme is: %f.\n", orderTipS);
    fprintf("The computed order of convergence for tip force in the cubic spline (local) scheme is: %f.\n", orderTipSL);
end