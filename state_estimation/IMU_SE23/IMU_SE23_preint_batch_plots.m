% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        SE2(3) with bias IMU Preintegration Simulated Data
% MODULE:       Plot
% DESCRIPTION:  Plots the true, noisy and MLG batch filtered data.
%               Additionally, plots the error in the state and
%               the innovation term with the 3 sigma bound for each
%               throughout the provided time.

%close all

font_size = 18;
line_size = 15;
line_width = 4;

T = gt_states.time(1:num_states);
num_rows = 5;
figure
subplot(num_xi/num_rows,num_rows,1);
plot (T,x_true(1,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(1,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,6);
plot (T,x_true(2,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(2,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,11);
plot (T,x_true(3,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(3,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,2);
plot (T,x_true(4,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(4,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,7);
plot (T,x_true(5,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(5,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,12);
plot (T,x_true(6,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(6,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,3);
plot (T,x_true(7,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(7,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,8);
plot (T,x_true(8,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(8,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,13);
plot (T,x_true(9,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(9,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,4);
plot (T,x_true(10,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(10,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,9);
plot (T,x_true(11,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(11,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,14);
plot (T,x_true(12,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(12,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,5);
plot (T,x_true(13,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(13,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_1^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,10);
plot (T,x_true(14,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(14,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_2^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,15);
plot (T,x_true(15,:),'k','Linewidth',line_width)
hold on
grid on
grid minor
plot (T,x_batch(15,:),'Color','[0.3010, 0.7450, 0.9330]','Linewidth',line_width)
%ylim([1.5,3.5])
%%leg = legend("True", "Batch Estimate",'fontsize',font_size);
%leg.NumColumns=2;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\xi_3^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

%%

figure
subplot(num_xi/num_rows,num_rows,1);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(1,:), fliplr(-sig_batch(1,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(1,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
leg = legend({'$\pm3\sigma$ Bound', 'Error in $\xi$'},'Interpreter','latex','fontsize',font_size);
leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,6);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(2,:), fliplr(-sig_batch(2,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(2,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,11);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(3,:), fliplr(-sig_batch(3,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(3,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^\phi$ (rad)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,2);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(4,:), fliplr(-sig_batch(4,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(4,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,7);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(5,:), fliplr(-sig_batch(5,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(5,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
%xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,12);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(6,:), fliplr(-sig_batch(6,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(6,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^v$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,3);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(7,:), fliplr(-sig_batch(7,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(7,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,8);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(8,:), fliplr(-sig_batch(8,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(8,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,13);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(9,:), fliplr(-sig_batch(9,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(9,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^r$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,4);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(10,:), fliplr(-sig_batch(10,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(10,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,9);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(11,:), fliplr(-sig_batch(11,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(11,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,14);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(12,:), fliplr(-sig_batch(12,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(12,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^{b_g}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,5);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(13,:), fliplr(-sig_batch(13,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(13,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_1^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,10);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(14,:), fliplr(-sig_batch(14,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(14,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_2^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(num_xi/num_rows,num_rows,15);
x2 = [T, fliplr(T)];
inBetween = [sig_batch(15,:), fliplr(-sig_batch(15,:))];
fill(x2, inBetween, [0, 0.75, 0.75]);
hold on
grid on
grid minor
plot (T,err_batch(15,:),'k','Linewidth',line_width)
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Error $\xi_3^{b_a}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

%%
figure
subplot(3,1,1);
plot(T,RMSE_total_LI,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,RMSE_total_RI,'r','Linewidth',line_width);
%ylim([1.5,3.5])
leg = legend("Left Invariant", "Right Invariant",'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Total RMSE','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(3,1,2);
plot(T,RMSE_yaw_LI*180/pi,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,RMSE_yaw_RI*180/pi,'r','Linewidth',line_width);
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Yaw RMSE (deg)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(3,1,3);
plot(T,RMSE_pos_LI,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,RMSE_pos_RI,'r','Linewidth',line_width);
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Position RMSE (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

%%
figure
subplot(3,1,1);
plot(T,NEES_total_LI,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,NEES_total_RI,'r','Linewidth',line_width);
%ylim([1.5,3.5])
leg = legend("Left Invariant", "Right Invariant",'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Total NEES','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(3,1,2);
plot(T,NEES_yaw_LI,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,NEES_yaw_RI,'r','Linewidth',line_width);
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Yaw NEES','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)

subplot(3,1,3);
plot(T,NEES_pos_LI,'b','Linewidth',line_width);
hold on
grid on
grid minor
plot(T,NEES_pos_RI,'r','Linewidth',line_width);
%ylim([1.5,3.5])
%leg = legend({'$\pm3\sigma$ Bound', '$x_{\textrm{true}} - \hat{x}$'},'Interpreter','latex','fontsize',font_size);
%leg.NumColumns=1;
%title('\textbf{Cart Position in Front of a Wall with Time}','fontsize',font_size,'Interpreter','latex') 
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('Position NEES','fontsize',font_size,'Interpreter','latex');
set(gca,'FontSize',font_size)
