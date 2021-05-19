% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Extended Kalman Filter
% MODULE:       Plot
% DESCRIPTION:  Plots the true, noisy and Kalman filtered data for the
%               position of a cart in front of a wall. Additionally, plots 
%               the error in the position and the innovation term with the
%               3 sigma bound for each throughout the provided time.

close all

if plot_ode
    x_true = x_ode;
    x_cor = x_ode_cor;
    v_true = v_ode;
    v_cor = v_ode_cor;
else
    x_true = x_sin;
    x_cor = x_sin_cor;
    v_true = v_sin;
    v_cor = v_sin_cor;
end

font_size = 32;
line_size = 15;
line_width = 4;

figure
subplot(3,1,1);
plot (T,x_true,'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_cor,'Color','[0.8500, 0.3250, 0.0980]','Linewidth',1)
plot (T,x_KF,'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Measured State", "EKF Estimate",'fontsize',font_size);
leg.NumColumns=3;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Position (m)','fontsize',font_size,'Interpreter','latex');

subplot(3,1,2);
x2 = [T, fliplr(T)];
inBetween = [r_norm, fliplr(-r_norm)];
fill(x2, inBetween, [0, 0.75, 0.75]);
ylim([-3,3])
hold on
grid on
grid minor
plot (T,x_true-x_KF,'k','Linewidth',line_width)
%title('\textbf{Error in Position with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error (m)','fontsize',font_size,'Interpreter','latex');
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;

subplot(3,1,3);
if ~(same_freq == 2)
    x2 = [T_inov, fliplr(T_inov)];
    inBetween = [inov_norm, fliplr(-inov_norm)];
    fill(x2, inBetween, [0, 0.75, 0.75]);
    %title('\textbf{Innovation Term with Time}','fontsize',font_size,'Interpreter','latex') 
    xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
    ylabel('Innovation Term','fontsize',font_size,'Interpreter','latex');
    hold on
    grid on
    grid minor
    plot (T_inov,rho_k,'k','Linewidth',2)
    ylim([-15,15])
    leg = legend({'$\pm3\sigma$ Bound', '$y_{\textrm{true}} - \check{y}$'},'Interpreter','latex','fontsize',font_size);
    leg.NumColumns=1;
end
plot_name = 'plots/';
if same_freq == 1
    plot_name = strcat(plot_name,'EKF_synchronous_');
else
    plot_name = strcat(plot_name,'EKF_asynchronous_');
end
if plot_ode == 1
    plot_name = strcat(plot_name,'mass_spring.eps');
else
    plot_name = strcat(plot_name,'sinusoidal.eps');
end
%exportfig(gcf,plot_name,'width',12,'Height',10,'fontmode','fixed','fontsize',font_size,'Color','cmyk','LineWidth',line_width);