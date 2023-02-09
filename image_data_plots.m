function t = image_data_plots(prefix)
    % run the simulations
    for name = {[prefix, '-1'], [prefix, '-2'], [prefix, '-3'], [prefix, '-4'], [prefix, '-5']}
        growth_anisotropic_image(['../cell-profiles/', name{1}, '.png'], name{1}, 12)
    end
    
    % load the data from the simulations
    load(['../cell-profiles/', prefix, '-1.mat']); 
    s1 = s; gam1 = gam; gam_s1 = gam_s; gam_theta1 = gam_theta; eps_s1 = eps_s; eps_theta1 = eps_theta; velocity1 = velocity; ks1 = ks; ktheta1 = ktheta; strainl1 = strainl; strainr1 = strainr;
    load(['../cell-profiles/', prefix, '-2.mat']);
    s2 = s; gam2 = gam; gam_s2 = gam_s; gam_theta2 = gam_theta; eps_s2 = eps_s; eps_theta2 = eps_theta; velocity2 = velocity; ks2 = ks; ktheta2 = ktheta; strainl2 = strainl; strainr2 = strainr;
    load(['../cell-profiles/', prefix, '-3.mat']);
    s3 = s; gam3 = gam; gam_s3 = gam_s; gam_theta3 = gam_theta; eps_s3 = eps_s; eps_theta3 = eps_theta; velocity3 = velocity; ks3 = ks; ktheta3 = ktheta; strainl3 = strainl; strainr3 = strainr;
    load(['../cell-profiles/', prefix, '-4.mat']);
    s4 = s; gam4 = gam; gam_s4 = gam_s; gam_theta4 = gam_theta; eps_s4 = eps_s; eps_theta4 = eps_theta; velocity4 = velocity; ks4 = ks; ktheta4 = ktheta; strainl4 = strainl; strainr4 = strainr;
    load(['../cell-profiles/', prefix, '-5.mat']);
    s5 = s; gam5 = gam; gam_s5 = gam_s; gam_theta5 = gam_theta; eps_s5 = eps_s; eps_theta5 = eps_theta; velocity5 = velocity; ks5 = ks; ktheta5 = ktheta; strainl5 = strainl; strainr5 = strainr;
    
    s_max = max([s1 s2 s3 s4 s5]);
    
    % make some plots
    % \gamma
    f = figure;
    plot(s1, gam1/max(gam1), s2, gam2/max(gam2), s3, gam3/max(gam3), s4, gam4/max(gam4), s5, gam5/max(gam5), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma');
    pbaspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_gamma_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_gamma_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    % \gamma_s
    f = figure;
    plot(s1, gam_s1/max([gam_s1 gam_theta1]), s2, gam_s2/max([gam_s2 gam_theta2]), s3, gam_s3/max([gam_s3 gam_theta3]), s4, gam_s4/max([gam_s4 gam_theta4]), s5, gam_s5/max([gam_s5 gam_theta5]), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma_s');
    pbaspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_gamma_s_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_gamma_s_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    % \gamma_theta
    f = figure;
    plot(s1, gam_theta1/max([gam_s1 gam_theta1]), s2, gam_theta2/max([gam_s2 gam_theta2]), s3, gam_theta3/max([gam_s3 gam_theta3]), s4, gam_theta4/max([gam_s4 gam_theta4]), s5, gam_theta5/max([gam_s5 gam_theta5]), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma_\theta');
    pbaspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_gamma_theta_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_gamma_theta_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    eps_max = max([eps_s1 eps_theta1 eps_s2 eps_theta2 eps_s3 eps_theta3 eps_s4 eps_theta4 eps_s5 eps_theta5]);
    eps_min = min([eps_s1 eps_theta1 eps_s2 eps_theta2 eps_s3 eps_theta3 eps_s4 eps_theta4 eps_s5 eps_theta5]);
    eps_diff = eps_max - eps_min;
    
    % \epsilon_s
    f = figure;
    plot(s1(1:end-1), eps_s1/max([eps_s1 eps_theta1]), s2(1:end-1), eps_s2/max([eps_s2 eps_theta2]), s3(1:end-1), eps_s3/max([eps_s3 eps_theta3]), s4(1:end-1), eps_s4/max([eps_s4 eps_theta4]), s5(1:end-1), eps_s5/max([eps_s5 eps_theta5]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_s$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([-0.1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_epsilon_s_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_epsilon_s_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    % \epsilon_\theta
    f = figure;
    plot(s1(1:end-1), eps_theta1/max([eps_s1 eps_theta1]), s2(1:end-1), eps_theta2/max([eps_s2 eps_theta2]), s3(1:end-1), eps_theta3/max([eps_s3 eps_theta3]), s4(1:end-1), eps_theta4/max([eps_s4 eps_theta4]), s5(1:end-1), eps_theta5/max([eps_s5 eps_theta5]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_\theta$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([-0.1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_epsilon_theta_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_epsilon_theta_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    curv_max = max([ks1 ktheta1 ks2 ktheta2 ks3 ktheta3 ks4 ktheta4 ks5 ktheta5]);
    curv_min = min([ks1 ktheta1 ks2 ktheta2 ks3 ktheta3 ks4 ktheta4 ks5 ktheta5]);
    curv_diff = curv_max - curv_min;
    
    % \kappa_s
    f = figure;
    plot(s1, fliplr(ks1), s2, fliplr(ks2), s3, fliplr(ks3), s4, fliplr(ks4), s5, fliplr(ks5), 'LineWidth', 2.0);
    xlabel('s'); ylabel('\kappa_s');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([curv_min-0.1*curv_diff curv_max+0.1*curv_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_kappa_s_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_kappa_s_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    % \kappa_\theta
    f = figure;
    plot(s1, fliplr(ktheta1), s2, fliplr(ktheta2), s3, fliplr(ktheta3), s4, fliplr(ktheta4), s5, fliplr(ktheta5), 'LineWidth', 2.0);
    xlabel('s'); ylabel('\kappa_\theta');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([curv_min-0.1*curv_diff curv_max+0.1*curv_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_kappa_theta_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_kappa_theta_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    lambda_max = max([strainl1 strainr1 strainl2 strainr2 strainl3 strainr3 strainl4 strainr4 strainl5 strainr5]);
    lambda_min = min([strainl1 strainr1 strainl2 strainr2 strainl3 strainr3 strainl4 strainr4 strainl5 strainr5]);
    lambda_diff = lambda_max - lambda_min;
    
    % \lambda_s
    f = figure;
    plot(s1, fliplr(strainl1(2:end)), s2, fliplr(strainl2(2:end)), s3, fliplr(strainl3(2:end)), s4, fliplr(strainl4(2:end)), s5, fliplr(strainl5(2:end)), 'LineWidth', 2.0);
    xlabel('s'); ylabel('\lambda_s');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_lambda_s_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_lambda_s_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
    
    % \lambda_\theta
    f = figure;
    plot(s1, fliplr(strainr1(2:end)), s2, fliplr(strainr2(2:end)), s3, fliplr(strainr3(2:end)), s4, fliplr(strainr4(2:end)), s5, fliplr(strainr5(2:end)), 'LineWidth', 2.0);
    xlabel('s'); ylabel('\lambda_\theta');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, ['media/anisotropic_', prefix, '_lambda_theta_aggregate.png']);
    f.PaperPosition = [1.333 3.31 5 5];
    filename = ['media/anisotropic_', prefix, '_lambda_theta_aggregate.png'];
    print(filename, '-dpng', '-r100');
    close all;
end