function t = aggregate_hyphoid_data_plots
    stiffness = 24;
%     a_vals = [3.4 4.1 4.3 4.55 4.9 5.6 7.0 8.8]; % some hyphoid a-values, including those for moss chloronema and caulonema (3.4, 4.1) and root hairs (6.8)
    a_vals = [3.4 4.1 4.8 5.6 6.8 8.8]; % some hyphoid a-values, including those for moss chloronema and caulonema (3.4, 4.1) and root hairs (6.8)
    
    % run the simulations
%     for a = a_vals
%         growth_anisotropic_hyphoid(a, stiffness);
%     end
    
    % load the data from the simulations
    load(['hyphoid_a_', num2str(a_vals(1)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s1 = s; gam1 = gam; gam_s1 = gam_s; gam_theta1 = gam_theta; eps_s1 = eps_s; eps_theta1 = eps_theta; eps_se1 = eps_se; eps_thetae1 = eps_thetae; extensibility1 = extensibility; velocity1 = velocity; ks1 = ksspline; ktheta1 = kthetaspline; strainl1 = strainl; strainr1 = strainr;
    load(['hyphoid_a_', num2str(a_vals(2)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s2 = s; gam2 = gam; gam_s2 = gam_s; gam_theta2 = gam_theta; eps_s2 = eps_s; eps_theta2 = eps_theta; eps_se2 = eps_se; eps_thetae2 = eps_thetae; extensibility2 = extensibility; velocity2 = velocity; ks2 = ksspline; ktheta2 = kthetaspline; strainl2 = strainl; strainr2 = strainr;
    load(['hyphoid_a_', num2str(a_vals(3)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s3 = s; gam3 = gam; gam_s3 = gam_s; gam_theta3 = gam_theta; eps_s3 = eps_s; eps_theta3 = eps_theta; eps_se3 = eps_se; eps_thetae3 = eps_thetae; extensibility3 = extensibility; velocity3 = velocity; ks3 = ksspline; ktheta3 = kthetaspline; strainl3 = strainl; strainr3 = strainr;
    load(['hyphoid_a_', num2str(a_vals(4)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s4 = s; gam4 = gam; gam_s4 = gam_s; gam_theta4 = gam_theta; eps_s4 = eps_s; eps_theta4 = eps_theta; eps_se4 = eps_se; eps_thetae4 = eps_thetae; extensibility4 = extensibility; velocity4 = velocity; ks4 = ksspline; ktheta4 = kthetaspline; strainl4 = strainl; strainr4 = strainr;
    load(['hyphoid_a_', num2str(a_vals(5)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s5 = s; gam5 = gam; gam_s5 = gam_s; gam_theta5 = gam_theta; eps_s5 = eps_s; eps_theta5 = eps_theta; eps_se5 = eps_se; eps_thetae5 = eps_thetae; extensibility5 = extensibility; velocity5 = velocity; ks5 = ksspline; ktheta5 = kthetaspline; strainl5 = strainl; strainr5 = strainr;
    load(['hyphoid_a_', num2str(a_vals(6)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
    s6 = s; gam6 = gam; gam_s6 = gam_s; gam_theta6 = gam_theta; eps_s6 = eps_s; eps_theta6 = eps_theta; eps_se6 = eps_se; eps_thetae6 = eps_thetae; extensibility6 = extensibility; velocity6 = velocity; ks6 = ksspline; ktheta6 = kthetaspline; strainl6 = strainl; strainr6 = strainr;
%     load(['hyphoid_a_', num2str(a_vals(7)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
%     s7 = s; gam7 = gam; gam_s7 = gam_s; gam_theta7 = gam_theta; eps_s7 = eps_s; eps_theta7 = eps_theta; velocity7 = velocity; ks7 = ksspline; ktheta7 = kthetaspline; strainl7 = strainl; strainr7 = strainr;
%     load(['hyphoid_a_', num2str(a_vals(8)), '_k_', num2str(stiffness), '_mu_', num2str(stiffness), '.mat']); 
%     s8 = s; gam8 = gam; gam_s8 = gam_s; gam_theta8 = gam_theta; eps_s8 = eps_s; eps_theta8 = eps_theta; velocity8 = velocity; ks8 = ksspline; ktheta8 = kthetaspline; strainl8 = strainl; strainr8 = strainr;
    
    s_max = max([s1 s2 s3 s4 s5 s6]);
    
    % make some plots
    % shapes
    f = figure;
    hold on; 
    for a = a_vals
        plot(pi / a * (-0.98:0.001:0.98) .* cot(pi * (-0.98:0.001:0.98)) - 1/a, (-0.98:0.001:0.98), 'LineWidth', 2.0); 
    end
    daspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
    xlim([-3.2 0]);
    daspect([1 1 1]);
    pbaspect([1 1 1]);
%     set(gca,'xticklabel',[],'yticklabel',[])
    xlabel('z');
    ylabel('r');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
%     exportgraphics(f, 'media/hyphoid-aggregate-profiles.png');
    print -dpng 'media/hyphoid-aggregate-profiles.png' -r100
    close all;
    
    % \gamma
    f = figure;
    plot(s1, gam1/max(gam1), s2, gam2/max(gam2), s3, gam3/max(gam3), s4, gam4/max(gam4), s5, gam5/max(gam5), s6, gam6/max(gam6), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_gamma_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_gamma_aggregate.png' -r100
    close all;
    
    % \gamma_s
    f = figure;
    plot(s1, gam_s1/max([gam_s1 gam_theta1]), s2, gam_s2/max([gam_s2 gam_theta2]), s3, gam_s3/max([gam_s3 gam_theta3]), s4, gam_s4/max([gam_s4 gam_theta4]), s5, gam_s5/max([gam_s5 gam_theta5]), s6, gam_s6/max([gam_s6 gam_theta6]), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma_s');
    pbaspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_gamma_s_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_gamma_s_aggregate.png' -r100
    close all;
    
    % \gamma_theta
    f = figure;
    plot(s1, gam_theta1/max([gam_s1 gam_theta1]), s2, gam_theta2/max([gam_s2 gam_theta2]), s3, gam_theta3/max([gam_s3 gam_theta3]), s4, gam_theta4/max([gam_s4 gam_theta4]), s5, gam_theta5/max([gam_s5 gam_theta5]), s6, gam_theta6/max([gam_s6 gam_theta6]), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\gamma_\theta');
    pbaspect([1 1 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_gamma_theta_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_gamma_theta_aggregate.png' -r100
    close all;
    
    eps_max = max([eps_s1 eps_theta1 eps_s2 eps_theta2 eps_s3 eps_theta3 eps_s4 eps_theta4 eps_s5 eps_theta5 eps_s6 eps_theta6]);
    eps_min = min([eps_s1 eps_theta1 eps_s2 eps_theta2 eps_s3 eps_theta3 eps_s4 eps_theta4 eps_s5 eps_theta5 eps_s6 eps_theta6]);
    eps_diff = eps_max - eps_min;
    
    % \epsilon_s
    f = figure;
    plot(s1(1:end-1), eps_s1/max([eps_s1 eps_theta1]), s2(1:end-1), eps_s2/max([eps_s2 eps_theta2]), s3(1:end-1), eps_s3/max([eps_s3 eps_theta3]), s4(1:end-1), eps_s4/max([eps_s4 eps_theta4]), s5(1:end-1), eps_s5/max([eps_s5 eps_theta5]), s6(1:end-1), eps_s6/max([eps_s6 eps_theta6]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_s$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([0 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_epsilon_s_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_epsilon_s_aggregate.png' -r100
    close all;
    
    % \epsilon_\theta
    f = figure;
    plot(s1(1:end-1), eps_theta1/max([eps_s1 eps_theta1]), s2(1:end-1), eps_theta2/max([eps_s2 eps_theta2]), s3(1:end-1), eps_theta3/max([eps_s3 eps_theta3]), s4(1:end-1), eps_theta4/max([eps_s4 eps_theta4]), s5(1:end-1), eps_theta5/max([eps_s5 eps_theta5]), s6(1:end-1), eps_theta6/max([eps_s6 eps_theta6]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_\theta$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([0 1]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_epsilon_theta_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_epsilon_theta_aggregate.png' -r100
    close all;
    
    % \epsilon_se
    f = figure;
    plot(s1(1:end-1), eps_se1/max([eps_se1 eps_thetae1]), s2(1:end-1), eps_se2/max([eps_se2 eps_thetae2]), s3(1:end-1), eps_se3/max([eps_se3 eps_thetae3]), s4(1:end-1), eps_se4/max([eps_se4 eps_thetae4]), s5(1:end-1), eps_se5/max([eps_se5 eps_thetae5]), s6(1:end-1), eps_se6/max([eps_se6 eps_thetae6]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_{se}$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([-0.4 1.2]);
    ax = gca;
    set(ax, 'fontsize', 14);
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_epsilon_se_aggregate.png' -r100
    close all;
    
    % \epsilon_\thetae
    f = figure;
    plot(s1(1:end-1), eps_thetae1/max([eps_se1 eps_thetae1]), s2(1:end-1), eps_thetae2/max([eps_se2 eps_thetae2]), s3(1:end-1), eps_thetae3/max([eps_se3 eps_thetae3]), s4(1:end-1), eps_thetae4/max([eps_se4 eps_thetae4]), s5(1:end-1), eps_thetae5/max([eps_se5 eps_thetae5]), s6(1:end-1), eps_thetae6/max([eps_se6 eps_thetae6]), 'LineWidth', 2.0);
    xlabel('s'); ylabel('$\dot{\epsilon}_{\theta e}$', 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([-0.4 1.2]);
    ax = gca;
    set(ax, 'fontsize', 14);
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_epsilon_thetae_aggregate.png' -r100
    close all;
    
    % extensibility
    f = figure;
    plot(s1, extensibility1(1:end-1)/max(extensibility1), s2, extensibility2(1:end-1)/max(extensibility2), s3, extensibility3(1:end-1)/max(extensibility3), s4, extensibility4(1:end-1)/max(extensibility4), s5, extensibility5(1:end-1)/max(extensibility5), s6, extensibility6(1:end-1)/max(extensibility6), 'LineWidth', 2.0);
    % xlim([0 max(s1)]); ylim([0 1]);
    xlabel('s'); ylabel('\Phi');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ax = gca;
    set(ax, 'fontsize', 14);
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_extensibility_aggregate.png' -r100
    close all;
    
    curv_max = max([ks1 ktheta1 ks2 ktheta2 ks3 ktheta3 ks4 ktheta4 ks5 ktheta5 ks6 ktheta6]);
    curv_min = min([ks1 ktheta1 ks2 ktheta2 ks3 ktheta3 ks4 ktheta4 ks5 ktheta5 ks6 ktheta6]);
    curv_diff = curv_max - curv_min;
    
    % \kappa_s
    f = figure;
    hold on;
    plot(s1, fliplr(ks1(2:end)), s2, fliplr(ks2(2:end)), s3, fliplr(ks3(2:end)), s4, fliplr(ks4(2:end)), s5, fliplr(ks5(2:end)), s6, fliplr(ks6(2:end)), 'LineWidth', 2.0);
    plot(s1(end-find(ks1 == max(ks1))+2), ks1(find(ks1 == max(ks1))), 'xk', ...
        s2(end-find(ks2 == max(ks2))+2), ks2(find(ks2 == max(ks2))), 'xk', ...
        s3(end-find(ks3 == max(ks3))+2), ks3(find(ks3 == max(ks3))), 'xk', ...
        s4(end-find(ks4 == max(ks4))+2), ks4(find(ks4 == max(ks4))), 'xk', ...
        s5(end-find(ks5 == max(ks5))+2), ks5(find(ks5 == max(ks5))), 'xk', ...
        s6(end-find(ks6 == max(ks6))+2), ks6(find(ks6 == max(ks6))), 'xk', 'LineWidth', 2.0);
    xlabel('s'); ylabel('\kappa_s');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([0 curv_max+0.1*curv_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_kappa_s_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_kappa_s_aggregate.png' -r100
    close all;
    
    % \kappa_\theta
    f = figure;
    plot(s1, fliplr(ktheta1(2:end)), s2, fliplr(ktheta2(2:end)), s3, fliplr(ktheta3(2:end)), s4, fliplr(ktheta4(2:end)), s5, fliplr(ktheta5(2:end)), s6, fliplr(ktheta6(2:end)), 'LineWidth', 2.0);
    xlabel('s'); ylabel('\kappa_\theta');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([0 curv_max+0.1*curv_diff]);
    ax = gca;
    set(ax, 'fontsize', 14);
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_kappa_theta_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_kappa_theta_aggregate.png' -r100
    close all;
    
    lambda_max = max([strainl1 strainr1 strainl2 strainr2 strainl3 strainr3 strainl4 strainr4 strainl5 strainr5 strainl6 strainr6]);
    lambda_min = min([strainl1 strainr1 strainl2 strainr2 strainl3 strainr3 strainl4 strainr4 strainl5 strainr5 strainl6 strainr6]);
    lambda_diff = lambda_max - lambda_min;
    
    % \lambda_s
    f = figure;
%     plot(s1, fliplr(strainl1(2:end)), s2, fliplr(strainl2(2:end)), s3, fliplr(strainl3(2:end)), s4, fliplr(strainl4(2:end)), s5, fliplr(strainl5(2:end)), s6, fliplr(strainl6(2:end)), 'LineWidth', 2.0);
%     xlabel('s'); ylabel('\lambda_s');
%     pbaspect([1 1 1]);
%     xlim([0 s_max]);
%     ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff]);
    plot(s1, fliplr(strainl1(2:end)) - 1, s2, fliplr(strainl2(2:end)) - 1, s3, fliplr(strainl3(2:end)) - 1, s4, fliplr(strainl4(2:end)) - 1, s5, fliplr(strainl5(2:end)) - 1, s6, fliplr(strainl6(2:end)) - 1, 'LineWidth', 2.0);
    xlabel('s'); ylabel('\lambda_s - 1');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff] - 1);
    ax = gca;
    set(ax, 'fontsize', 14);
    ax.YAxis.Exponent = -2;
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_lambda_s_aggregate.png');
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_lambda_s_aggregate.png' -r100
    close all;
    
    % \lambda_\theta
    f = figure;
%     plot(s1, fliplr(strainr1(2:end)), s2, fliplr(strainr2(2:end)), s3, fliplr(strainr3(2:end)), s4, fliplr(strainr4(2:end)), s5, fliplr(strainr5(2:end)), s6, fliplr(strainr6(2:end)), 'LineWidth', 2.0);
%     xlabel('s'); ylabel('\lambda_\theta');
%     pbaspect([1 1 1]);
%     xlim([0 s_max]);
%     ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff]);
    plot(s1, fliplr(strainr1(2:end)) - 1, s2, fliplr(strainr2(2:end)) - 1, s3, fliplr(strainr3(2:end)) - 1, s4, fliplr(strainr4(2:end)) - 1, s5, fliplr(strainr5(2:end)) - 1, s6, fliplr(strainr6(2:end)) - 1, 'LineWidth', 2.0);
    xlabel('s'); ylabel('\lambda_\theta - 1');
    pbaspect([1 1 1]);
    xlim([0 s_max]);
    ylim([lambda_min-0.1*lambda_diff lambda_max+0.1*lambda_diff] - 1);
    ax = gca;
    set(ax, 'fontsize', 14);
    ax.YAxis.Exponent = -2;
%     exportgraphics(gcf, 'media/anisotropic_hyphoid_lambda_theta_aggregate.png');
%     f.PaperPosition = [1.333 3.31 5 5];
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    set(gca,'Box','on');
    pause(0.5);
    print -dpng 'media/anisotropic_hyphoid_lambda_theta_aggregate.png' -r100
    close all;
end