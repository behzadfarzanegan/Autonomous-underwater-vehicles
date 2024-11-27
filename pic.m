
clc
clear all
close all
%%%%%%%%%%%%%%%
figure(1);
load('initial_control.mat')
plot(eta(1,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(eta(1,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(eta(1,:),'-.b','LineWidth',2)
hold on

plot(eta_d(1,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Surge ({x}) [m]','FontWeight','b','FontSize', 12 );
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10,'location', 'northwest')
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.53 .25 .33 .25]);
plot(ax1,eta(1,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax1,eta_d(1,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 10000])

set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
load('initial_control.mat')
plot(eta(2,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(eta(2,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(eta(2,:),'-.b','LineWidth',2)
hold on

plot(eta_d(2,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Sway (y) [m]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10,'NumColumns',2 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );

%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
load('initial_control.mat')
plot(eta(3,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(eta(3,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(eta(3,:),'-.b','LineWidth',2)
hold on

plot(eta_d(3,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Heave (z) [m]','FontWeight','b','FontSize', 12) ;
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10 ,'Location','north','NumColumns',2)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.2 .26  .33 .25]);
plot(ax1,eta(3,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax1,eta_d(3,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 3000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

ax2 = axes('Position',[.6 .46  .28 .25]);
plot(ax2,eta(3,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax2,eta_d(3,:),'--r','LineWidth',2);
hold on
load('without_LL.mat')
plot(ax2,eta(3,:),':g','LineWidth',2);
ax2.XAxis.Exponent=0;
xlim([45000 49999])
set( gca, 'FontWeight', 'b','FontSize', 10 );
set(gca,'XTick',[], 'YTick', [])

% ax1=gca;
% ax1.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
load('initial_control.mat')
plot(eta(4,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(eta(4,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(eta(4,:),'-.b','LineWidth',2)
hold on

plot(eta_d(4,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Pitch (\theta) [rad]','FontWeight','b','FontSize', 12 );
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10,'Location','southeast','NumColumns',2 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax2 = axes('Position',[.55 .46  .3 .25]);
plot(ax2,eta(4,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax2,eta_d(4,:),'--r','LineWidth',2);
hold on
load('without_LL.mat')
plot(ax2,eta(4,:),':g','LineWidth',2);
ax2.XAxis.Exponent=0;
xlim([30000 50000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
load('initial_control.mat')
plot(eta(5,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(eta(5,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(eta(5,:),'-.b','LineWidth',2)
hold on

plot(eta_d(5,:),'--r','LineWidth',2)

xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Yaw (\psi) [rad]','FontWeight','b','FontSize', 12 );
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10,'NumColumns',2 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%

figure(6);
load('initial_control.mat')
plot3(eta(1,:),eta(2,:),eta(3,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot3(eta(1,:),eta(2,:),eta(3,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot3(eta(1,:),eta(2,:),eta(3,:),'-.b','LineWidth',2)

hold on

plot3(eta_d(1,:),eta_d(2,:),eta_d(3,:),'--r','LineWidth',2)
view(-12,17)

xlabel('x [m]','FontWeight','b','FontSize', 14);
ylabel('y [m]','FontWeight','b','FontSize', 14);
zlabel('z [m]','FontWeight','b','FontSize', 14);
legend('PD', 'Without LL','With LL', 'Reference','FontWeight','b','FontSize', 10,'location','northeast')
set( gca, 'FontWeight', 'b','FontSize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(7);
load('initial_control.mat')
plot(x2(1,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(x2(1,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(x2(1,:),'-.b','LineWidth',2)
hold on

plot(V1(1,:),'--r','LineWidth',2)

xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Surge velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL',  'Virtual control','FontWeight','b','FontSize', 10,'NumColumns',2)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12);
ax=gca;
ax.XAxis.Exponent=2;

ax2 = axes('Position',[.55 .36  .3 .25]);
plot(ax2,x2(1,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax2,V1(1,:),'--r','LineWidth',2);
hold on
load('without_LL.mat')
plot(ax2,x2(1,:),':g','LineWidth',2);

hold on
load('initial_control.mat')
plot(ax2,x2(1,:),'k','LineWidth',2);

ax2.XAxis.Exponent=0;
xlim([30000 50000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8);
load('initial_control.mat')
plot(x2(2,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(x2(2,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(x2(2,:),'-.b','LineWidth',2)
hold on

plot(V1(2,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Sway velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL',  'Virtual control','FontWeight','b','FontSize', 10,'NumColumns',2)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

% ax2 = axes('Position',[.2 .6  .3 .25]);
% plot(ax2,x2(2,:),'-.b','LineWidth',2);
% grid off;box on
% hold on
% plot(ax2,V1(2,:),'--r','LineWidth',2);
% hold on
% load('without_LL.mat')
% plot(ax2,x2(2,:),':g','LineWidth',2);
% 
% hold on
% load('initial_control.mat')
% plot(ax2,x2(2,:),'k','LineWidth',2);
% 
% ax2.XAxis.Exponent=0;
% xlim([45000 50000])
% set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%


figure(9);
load('initial_control.mat')
plot(x2(3,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(x2(3,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(x2(3,:),'-.b','LineWidth',2)
hold on

plot(V1(3,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Heave velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL',  'Virtual control','FontWeight','b','FontSize', 9,'NumColumns',2)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.3 .26  .4 .3]);
plot(ax1,x2(3,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax1,V1(3,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 4000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);
load('initial_control.mat')
plot(x2(4,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(x2(4,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(x2(4,:),'-.b','LineWidth',2)
hold on

plot(V1(4,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Pitch velocity [rad.s^{-1}]','FontWeight','b','FontSize', 12 );
legend('PD', 'Without LL','With LL',  'Virtual control','FontWeight','b','FontSize', 10,'NumColumns',2 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.24 .44  .38 .28]);
plot(ax1,x2(4,:),'-.b','LineWidth',2);
grid off;box on
hold on
plot(ax1,V1(4,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 4000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);
load('initial_control.mat')
plot(x2(5,:),'k','LineWidth',2)
hold on

load('without_LL.mat')
plot(x2(5,:),':g','LineWidth',2)
hold on

load('with_LL.mat')
plot(x2(5,:),'-.b','LineWidth',2)
hold on


plot(V1(5,:),'--r','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('Yaw velocity [rad.s^{-1}]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL',  'Virtual control','FontWeight','b','FontSize', 10,'NumColumns',2)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12);
subplot(3,1,1);
load('initial_control.mat')
plot(u(1,:),'k','LineWidth',2)
hold on 

load('without_LL.mat')
plot(u(1,:),':g','LineWidth',2)
hold on 

load('with_LL.mat')
plot(u(1,:),'-.b','LineWidth',2)
hold on 

xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_1 [N]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL','FontWeight','b','FontSize', 12)
xlim([0 NF])
ylim([-100 1000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.2 .81  .4 .1]);
plot(ax1,u(1,:),'-.b','LineWidth',2);
grid off;box on
ax1.XAxis.Exponent=0;
xlim([5000 50000])
set( gca, 'FontWeight', 'b','FontSize', 12 );

subplot(3,1,2);
load('initial_control.mat')
plot(u(2,:),'k','LineWidth',2)
hold on 

load('without_LL.mat')
plot(u(2,:),':g','LineWidth',2)
hold on 

load('with_LL.mat')
plot(u(2,:),'-.b','LineWidth',2)
hold on 
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_2 [N.m]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL','FontWeight','b','FontSize', 10)
xlim([0 NF])
ylim([-1000 800])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

subplot(3,1,3);
load('initial_control.mat')
plot(u(3,:),'k','LineWidth',2)
hold on 

load('without_LL.mat')
plot(u(3,:),':g','LineWidth',2)
hold on 

load('with_LL.mat')
plot(u(3,:),'-.b','LineWidth',2)
hold on 

xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_3 [N.m]','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL','FontWeight','b','FontSize', 12)
xlim([0 NF])
ylim([-100 1000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
ax1 = axes('Position',[.2 .21  .4 .1]);
plot(ax1,u(3,:),'-.b','LineWidth',2);
grid off;box on
ax1.XAxis.Exponent=0;
xlim([5000 50000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
%%%%%%%%%%%%%%%%%%%%%

figure(13);
% hold on 
load('initial_control')
plot(cost_com(:),'k','LineWidth',2)
hold on
load('without_LL')
plot(cost_com(:),':g','LineWidth',2)
hold on
load('with_LL')
plot(cost_com(:),'-.b','LineWidth',2)


xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('Cost','FontWeight','b','FontSize', 12);
legend('PD', 'Without LL','With LL','FontWeight','b','FontSize', 10)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.23 .41  .38 .28]);

plot(ax1,cost_com(:),'-.b','LineWidth',2);
hold on
load('without_LL')
plot(cost_com(:),':g','LineWidth',2)
grid off;box on
ax1.XAxis.Exponent=0;
xlim([5000 50000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
%%%%%%%%%%%%%%%%%%%%




for i=1:13
    h=figure(i);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4 3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 4 3]);
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,"fig"+num2str(i),'-dpdf','-r300')
     print(h,"fig"+num2str(i),'-dpng','-r300')
end

h=figure(6);
exportgraphics(h,'fig6.pdf' ,'Resolution',300)
% %%%%%%%%%%%%%%%%%%
    h=figure(12);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4 3]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 4 3]);
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,"fig12",'-dpdf','-r300')
    print(h,"fig12",'-dpng','-r300')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
