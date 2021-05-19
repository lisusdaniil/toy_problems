% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       Plot
% DESCRIPTION:  Plots the true, noisy and MLG batch filtered data for the
%               position and velocity of a cart in front of a wall.
%               Additionally, plots the error in the position, velocity and
%               the innovation term with the 3 sigma bound for each
%               throughout the provided time.

close all

font_size = 18;
line_size = 15;
line_width = 4;

figure
subplot(6,1,1);
plot (T,x_true(1,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(1,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,2);
plot (T,x_true(2,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(2,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,3);
plot (T,x_true(3,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(3,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,4);
plot (T,x_true(4,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(4,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,5);
plot (T,x_true(5,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(5,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,6);
plot (T,x_true(6,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(6,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend("True", "Batch Estimate",'fontsize',font_size);
leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)


figure
subplot(6,1,1);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(1,:), fliplr(-sig_batch(1,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(1,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,2);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(2,:), fliplr(-sig_batch(2,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(2,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,3);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(3,:), fliplr(-sig_batch(3,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(3,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,4);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(4,:), fliplr(-sig_batch(4,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(4,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,5);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(5,:), fliplr(-sig_batch(5,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(5,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(6,1,6);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(6,:), fliplr(-sig_batch(6,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(6,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

