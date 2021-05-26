% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear Batch Estimation
% MODULE:       Plot
% DESCRIPTION:  Plots the true, noisy and nonlinear batch filtered data for the
%               position and velocity of a cart in front of a wall.
%               Additionally, plots the error in the position, velocity and
%               the innovation term with the 3 sigma bound for each
%               throughout the provided time.

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
subplot(2,1,1);
plot (T,x_true,'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_cor,'Color','[0.8500, 0.3250, 0.0980]','Linewidth',1)
plot (T,x_batch,'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Measurement", "SWF Estimate",'fontsize',font_size);
leg.NumColumns=3;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Position (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(2,1,2);
x2 = [T, fliplr(T)];
inBetween = [3*sqrt(diag(Sigma)'), fliplr(-3*sqrt(diag(Sigma)'))];
fill(x2, inBetween, [0, 0.75, 0.75]);
%ylim([-3,3])
hold on
grid on
grid minor
plot (T,x_true-x_batch,'k','Linewidth',line_width)
%title('\textbf{Error in Position with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error (m)','fontsize',font_size,'Interpreter','latex');
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
set(gca,'FontSize',font_size)