
clear all;
clc;
 close all
rng(1)

eta(1,1)=1.2;
eta(1,2)=eta(1,1);

eta(2,1)=0;
eta(2,2)=eta(2,1);

eta(3,1)=42;
eta(3,2)=eta(3,1);

eta(4,1)=0;
eta(4,2)=eta(4,1);

eta(5,1)=0;
eta(5,2)=eta(5,1);
%%%%%%%%%%%%%%%%%%%
etak(1,1)=10.2;
etak(1,2)=etak(1,1);

etak(2,1)=0.2;
etak(2,2)=etak(2,1);

etak(3,1)=42;
etak(3,2)=etak(3,1);

etak(4,1)=0.5;
etak(4,2)=etak(4,1);

etak(5,1)=5;
etak(5,2)=etak(5,1);


%%%%%%%%%%%%%%%%%%%%%%%%%
eta_1d(1)=0;
eta_2d(1)=0;
eta_3d(1)=0;
eta_4d(1)=0;
eta_5d(1)=1;


eta_1d(2)=0;
eta_2d(2)=0;
eta_3d(2)=0;
eta_4d(2)=0;
eta_5d(2)=1;

eta_d(:,1)=[1 0 40 0 0]';
eta_d(:,2)=[1 0 40 0 0]';

etad(:,1)=[1 0 40 0 0]';
etad(:,2)=[1 0 40 0 0]';
%%%%%%%%%%%%%%%%%%%%%%
SS(:,1)=[0 0 0 0 0]';
SS(:,2)=[0 0 0 0 0]';

V1(:,1)=[0 0 0 0 0]';
V1(:,2)=[0 0 0 0 0]';
%%%%%%%%%%%%%%%%%%%%%
% uu(1)=0;v(1)=0;r(1)=0;
% uu(2)=0; v(2)=0; r(2)=0;
% nue(:,1)=[0 0 0]';
% nue(:,2)=[0 0 0]';

x(:,1) = [eta(1,1)-eta_1d(1),eta(2,1)-eta_2d(1),eta(3,1)-eta_3d(1),eta(4,1)-eta_4d(1),eta(5,1)-eta_5d(1)]'; %error
x(:,2) = [eta(1,1)-eta_1d(1),eta(2,1)-eta_2d(1),eta(3,1)-eta_3d(1),eta(4,1)-eta_4d(1),eta(5,1)-eta_5d(1)]'; 




Neuron_Num_a1=3;
Neuron_Num_c1=4;
Vc1=[rand(Neuron_Num_c1,2)];%no*1
Wc1=rand(1,Neuron_Num_c1)*.01;
 
Wc1(1,:)=Wc1;
Wc1(2,:)=Wc1(1,:);
Waa=0.01*rand(1,Neuron_Num_a1);
Wa1(1,:)=Waa;
Wa1(2,:)=Wa1(1,:);
r1=1;q1=100;ro1=15.1;kappa1=10.1;a1=6;

Neuron_Num_a2=3;
Neuron_Num_c2=4;
Vc2=[rand(Neuron_Num_c2,2)];%no*1
Wcc=rand(1,Neuron_Num_c2)*.01;
Wc2(1,:)=Wcc;
Wc2(2,:)=Wc2(1,:);

Wa2(1,:)=0.01*rand(1,Neuron_Num_a2);
Wa2(2,:)=Wa2(1,:);
r2=.1;q2=100;ro2=1.2;kappa2=10.1;a2=1;


Neuron_Num_a3=3;
Neuron_Num_c3=4;
Vc3=[rand(Neuron_Num_c3,2)];%no*1
Wc3(1,:)=rand(1,Neuron_Num_c3)*.01;
Wc3(2,:)=Wc3(1,:);
Wa3(1,:)=0.01*rand(1,Neuron_Num_a3);
Wa3(2,:)=Wa3(1,:);
r3=.1;q3=100;ro3=2.2;kappa3=10.7;a3=2;

Neuron_Num_a4=3;
Neuron_Num_c4=4;
Vc4=[rand(Neuron_Num_c4,2)];%no*1
Wc4(1,:)=rand(1,Neuron_Num_c4)*.01;
Wc4(2,:)=Wc4(1,:);
Wa4(1,:)=0.01*rand(1,Neuron_Num_a4);
Wa4(2,:)=Wa4(1,:);
r4=.1;q4=100;ro4=2.2;kappa4=10.7;a4=2;

Neuron_Num_a5=3;
Neuron_Num_c5=4;
Vc5=[rand(Neuron_Num_c5,2)];%no*1
Wc5(1,:)=rand(1,Neuron_Num_c5)*.01;
Wc5(2,:)=Wc5(1,:);
Wa5(1,:)=0.01*rand(1,Neuron_Num_a5);
Wa5(2,:)=Wa5(1,:);
r5=.1;q5=100;ro5=2.2;kappa5=10.7;a5=2;



gamma=0.99;
u1(1)= Wa1(1,:)*[x(1,1) x(1,1)^2 x(1,1)^3]'-1*x(1,1);
u1(2)= u1(1);

u2(1)= Wa2(1,:)*[x(2,1) x(2,1)^2 x(2,1)^3]'-1*x(2,1);
u2(2)= u2(1);

u3(1)= Wa3(1,:)*[x(3,1) x(3,1)^2 x(3,1)^3]'-1*x(3,1);
u3(2)= u3(1);

u4(1)= Wa4(1,:)*[x(4,1) x(4,1)^2 x(4,1)^3]'-1*x(4,1);
u4(2)= u4(1);

u5(1)= Wa5(1,:)*[x(5,1) x(5,1)^2 x(5,1)^3]'-1*x(5,1);
u5(2)= u5(1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2(:,1) = [0,0,0,1,1]';
x2(:,2) = [0,0,0,1,1]';
r(:,1) = [10,0,0,0,0]';
r(:,2) = [10,0,0,0,0]';
e(:,1)=x2(:,1)-r(:,1);
e(:,2)=x2(:,2)-r(:,2);
X(:,1) = [e(:,1);r(:,1);1];
X(:,2) = [e(:,2);r(:,2);1];
AugX(:,1) = [e(:,1);r(:,1)];
AugX(:,2) = [e(:,2);r(:,2)];


n=length(x);%states
m=3;%inputs
Q = 100000*eye(n);
R = .1*eye(m);

EWC=false;
Barrier = true;
BF=0;
cost_com(1)=0;
SS1(1)=0;
Fisher1=0;
Fisher2=0;
svMaxing = true;
expReplay = true;
%Z= zeros(2*n,2);p=0; %samples in replay buffer
buffSize = 36; %size of the buffer
dphi=[];

% Creat mesh domain
% aj =0.05;

H=[0 1 0 0 0 0]';
%au = 0.001;


%Neuron_Num_c = sum(1:length(AugX));
Neuron_Num_c = 15;
Neuron_Num_a =5;
%second layer weight update.
S = 300;

Q11=[];
Q22=[];
Q33=[];
Q44=[];

%%
%   bpTrain
% load('WEIHGT.mat')
%  InitializeWeights;

%    load('WcVc.mat')

Va = 2*randn(n,Neuron_Num_a);
Wa = -100*rand(Neuron_Num_a,m);
% load('WaVa.mat')
v_actor(:,:,1)=Va;
W_actor(:,:,1)=Wa;
v_actor(:,:,2)=v_actor(:,:,1);
W_actor(:,:,2)=W_actor(:,:,1);
%%%%%%%%%%%%
Wa1 = [-1000 +1000 1000 -1000 +1000]';
W_actor1(:,1)=Wa1;
W_actor1(:,2)=W_actor1(:,1);

Wa2 =  [-1000 1000 3000 -5000 -1000]';
W_actor2(:,1)=Wa2;
W_actor2(:,2)=W_actor2(:,1);

Wa3 = [-100 -100 1000 -100 -4000]';
W_actor3(:,1)=Wa3;
W_actor3(:,2)=W_actor3(:,1);

au1 = 0.02*0;
au2 = 0.02*0;
au3 = 0.02*0;
%%%%%%%%%%%%%%%
 % u(:,k+1) =[-15000*e(1,k)+1000*e(3,k) -5000*e(4,k) -4000*e(5,k)+1000*e(3,k)]';

%%%%%%%%%%%%%%%%%%
Wc = 100000*rand(Neuron_Num_c,1);
 %load('Wc.mat')
W_critic(:,:,1)=Wc;
W_critic(:,:,2)=W_critic(:,:,1);

aj =.2;
au = 0.002*1;
L=10;PE_C=1*0;M1=.5;

T=0.01;
NF=50000;

%       u(1,2) = W_actor1(:,2)'*[e(1,2) e(2,2) e(3,2) e(4,2) e(5,2)]';
%       u(2,2) = W_actor2(:,2)'*[e(1,2) e(2,2) e(3,2) e(4,2) e(5,2)]';
%       u(3,2) = W_actor3(:,2)'*[e(1,2) e(2,2) e(3,2) e(4,2) e(5,2)]';


ad_u = W_actor(:,:,2)'*logsig(v_actor(:,:,2)'*e(:,2));
% ad_u=[0.7;0.7];
u(:,1) = ad_u;
u(:,2) = ad_u;

Jhatsum(1:2)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:NF

[eta_d(:,k+1),r(:,k+1)] = Desired(eta_d(:,k), r(:,k),k);


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    R1=[cos(eta(5,k))*cos(eta(4,k)) -sin(eta(5,k)) sin(eta(4,k))*cos(eta(5,k))
        sin(eta(5,k))*cos(eta(4,k)) cos(eta(5,k)) sin(eta(4,k))*sin(eta(5,k))
        -sin(eta(4,k))              0             cos(eta(4,k))];

    R2=[1 0
        0 1/cos(eta(4,k))];

        %%%%%%%%%%%%%%%%%
     V1(1:3,k+1)=inv(R1)*(eta_d(1:3,k+1)-eta(1:3,k)-.999*(eta_d(1:3,k)-eta(1:3,k))+.01*([0  0 eta_d(3,k)]'-[0 0 eta(3,k)]'))/T;
     V1(4:5,k+1)=inv(R2)*(eta_d(4:5,k+1)-eta(4:5,k)-.99*(eta_d(4:5,k)-eta(4:5,k)))/T;
% e1=eta(:,k)-eta_d(:,k);
%      SS(1:3,k+1)=eta(1:3,k)+T*R1*V1(1:3,k)-eta_d(1:3,k+1);
%      SS(4:5,k+1)=eta(4:5,k)+T*R2*V1(4:5,k)-eta_d(4:5,k+1);
    
%     eta(1:3,k+1)=SS(1:3,k+1)+eta_d(1:3,k+1);
%     eta(4:5,k+1)=SS(4:5,k+1)+eta_d(4:5,k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%
  Yk_1 = e(:,k-1)'*Q*e(:,k-1)+u(:,k-1)'*R*u(:,k-1);%+1000*BF1 +10000*BF2;
    delta_x = gamma*Critic_NL_gamma_bah(e(:,k))-Critic_NL_gamma_bah(e(:,k-1));
    EJK = Yk_1+W_critic(:,:,k)'*delta_x;
    YK = e(:,k-1)'*Q*e(:,k-1)+u(:,k-1)'*R*u(:,k-1);%+1000*BF1 +10000*BF2;

    YY=e(:,k-1)'*Q*e(:,k-1)+u(:,k-1)'*R*u(:,k-1);
    cost_com(k)=cost_com(k-1)+YY;
    
%     Jhat(k) = W_critic(:,k)'*Critic_NL_gamma_bah(Inner(:,k));
%     Jhatsum(k+1)=Jhatsum(k)+Jhat(k);

    

    %%
     if k<20000
        Hessian=delta_x*delta_x';
        Fisher1=Fisher1+diag(diag(Hessian));
     elseif  k < 30000
        Hessian=delta_x*delta_x';
        Fisher2=Fisher2+diag(diag(Hessian));
      end
    temp_1 = W_critic(:,:,k);
  

    for LL=1:L
        delta_xK = gamma*Critic_NL_gamma_bah(e(:,k))-Critic_NL_gamma_bah(e(:,k-1));
        XK = delta_xK;
%                 if LL ==1
        if 20000 <= k && k < 30000

            temp_1 =temp_1-aj*(delta_xK*EJK)/(delta_xK'*delta_xK+1)-.000*aj*Fisher1*(temp_1-W_critic(:,:,19999))*EWC;
        elseif 30000<=k
            temp_1 =temp_1-aj*(delta_xK*EJK)/(delta_xK'*delta_xK+1)-0*aj*(temp_1-W_critic(:,:,19999))*EWC;
        else
            temp_1 =temp_1-aj*(delta_xK*EJK)/(delta_xK'*delta_xK+1);
        end

    end
    W_critic(:,:,k+1)=temp_1;


    [e(:,k+1),Gx] = NLinear_sys_NL_gamma_bah(x2(:,k),u(:,k),V1(:,k+1),k,eta(:,k));
    x2(:,k+1)=e(:,k+1)+V1(:,k+1);
    X(:,k+1) = [e(:,k+1);V1(:,k+1);1];
    AugX(:,k+1) = [e(:,k+1);V1(:,k+1)];

    uhat(:,k+1)= -gamma*0.5*inv(R)*Gx'*Dsigdx_NL_gamma_bah(e(:,k))'*W_critic(:,:,k);

    u_tilda1 = u(1,k)-uhat(1,k+1);
    u_tilda2 = u(2,k)-uhat(2,k+1);
    u_tilda3 = u(3,k)-uhat(3,k+1);

%        W_actor1(:,k+1) = W_actor1(:,k)-au1*W_actor1(:,k)*u_tilda1;
%        W_actor2(:,k+1) = W_actor2(:,k)-au2*W_actor2(:,k)*u_tilda2;
%        W_actor3(:,k+1) = W_actor3(:,k)-au3*W_actor3(:,k)*u_tilda3;
      
       phi=[e(1,k+1) e(2,k+1) e(3,k+1) e(4,k+1) e(5,k+1)]';
PHI=phi/(1+phi'*phi);


if 20000>k
        W_actor1(:,k+1) = W_actor1(:,k)-au1*PHI*u_tilda1;
       W_actor2(:,k+1) = W_actor2(:,k)-au2*PHI*u_tilda2;
       W_actor3(:,k+1) = W_actor3(:,k)-au3*PHI*u_tilda3;

   elseif 30000<=k
        W_actor1(:,k+1) = W_actor1(:,k)-au1*PHI*u_tilda1-20*au1*(W_actor1(:,k)-W_actor1(:,29999))*EWC;
       W_actor2(:,k+1) = W_actor2(:,k)-au2*PHI*u_tilda2-20*au2*(W_actor2(:,k)-W_actor2(:,29999))*EWC;
       W_actor3(:,k+1) = W_actor3(:,k)-au3*PHI*u_tilda3-20*au3*(W_actor3(:,k)-W_actor3(:,29999))*EWC;

else
        W_actor1(:,k+1) = W_actor1(:,k)-au1*PHI*u_tilda1-20*au1*(W_actor1(:,k)-W_actor1(:,19999))*EWC;
       W_actor2(:,k+1) = W_actor2(:,k)-au2*PHI*u_tilda2-20*au2*(W_actor2(:,k)-W_actor2(:,19999))*EWC;
       W_actor3(:,k+1) = W_actor3(:,k)-au3*PHI*u_tilda3-20*au3*(W_actor3(:,k)-W_actor3(:,19999))*EWC;

end 
      u(1,k+1) = W_actor1(:,k+1)'*[e(1,k+1) e(2,k+1) e(3,k+1) e(4,k+1) e(5,k+1)]';
      u(2,k+1) = W_actor2(:,k+1)'*[e(1,k+1) e(2,k+1) e(3,k+1) e(4,k+1) e(5,k+1)]';
      u(3,k+1) = W_actor3(:,k+1)'*[e(1,k+1) e(2,k+1) e(3,k+1) e(4,k+1) e(5,k+1)]';
 

%          u(:,k+1) =[-1500*e(1,k) -500*e(4,k)-500*e(3,k) -4000*e(5,k)]';
%        u(:,k+1) =[-15000*e(1,k)+1000*e(3,k) -5000*e(4,k) -4000*e(5,k)+1000*e(3,k)]';
%     tau=[5*d11-(m22*v*r-m33*w*q) -(m33-m11)*u*w-d55*q-69.42*sin(theta)+m55*(-theta +0.2-(thetad-theta)/T) -(m11-m22)*u*v]';
  
    U_tilda(:,k)=[u_tilda1 u_tilda2 u_tilda3]';
    Err(k)=EJK;
k
eta(:,k+1) = kinematic(eta(:,k),x2(:,k));
etad(:,k+1) = kinematic(etad(:,k),V1(:,k));
%%%%%%%;
end

figure(1);
hold on
plot(eta(1,:),'b','LineWidth',2)
hold on
plot(eta_d(1,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Surge ({x}) [m]','FontWeight','b','FontSize', 12 );
legend('Actual', 'Reference','FontWeight','b','FontSize', 12,'location', 'northwest')
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.53 .25 .33 .25]);
plot(ax1,eta(1,:),'b','LineWidth',2);
grid off;box on
hold on
plot(ax1,eta_d(1,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 10000])

set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
hold on
plot(eta(2,:),'b','LineWidth',2)
hold on
plot(eta_d(2,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Sway (y) [m]','FontWeight','b','FontSize', 12);
legend('Actual', 'Reference','FontWeight','b','FontSize', 12 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
hold on
plot(eta(3,:),'b','LineWidth',2)
hold on
plot(eta_d(3,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Heave (z) [m]','FontWeight','b','FontSize', 12) ;
legend('Actual', 'Reference','FontWeight','b','FontSize', 12 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.2 .26  .33 .25]);
plot(ax1,eta(3,:),'b','LineWidth',2);
grid off;box on
hold on
plot(ax1,eta_d(3,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 3000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

% ax1=gca;
% ax1.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
hold on
plot(eta(4,:),'b','LineWidth',2)
hold on
plot(eta_d(4,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Pitch (\theta) [rad]','FontWeight','b','FontSize', 12 );
legend('Actual', 'Reference','FontWeight','b','FontSize', 12 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
hold on
plot(eta(5,:),'b','LineWidth',2)
hold on
plot(eta_d(5,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Yaw (\psi) [rad]','FontWeight','b','FontSize', 12 );
legend('Actual', 'Reference','FontWeight','b','FontSize', 12 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%

figure(6);
hold on
plot3(eta(1,:),eta(2,:),eta(3,:),'b','LineWidth',2)
hold on
plot3(eta_d(1,:),eta_d(2,:),eta_d(3,:),'--r','LineWidth',2)
view(-12,17)

xlabel('x [m]','FontWeight','b','FontSize', 12);
ylabel('y [m]','FontWeight','b','FontSize', 12);
zlabel('z [m]','FontWeight','b','FontSize', 12);
legend('Actual', 'Reference','FontWeight','b','FontSize', 12,'location','northeast')
set( gca, 'FontWeight', 'b','FontSize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(7);
hold on
plot(x2(1,:),'b','LineWidth',2)
hold on
plot(V1(1,:),'--r','LineWidth',2)

xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Surge velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('Actual', 'Virtual control','FontWeight','b','FontSize', 12)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12);
ax=gca;
ax.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8);
hold on 
plot(x2(2,:),'b','LineWidth',2)
hold on
plot(V1(2,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Sway velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('Actual', 'Virtual control','FontWeight','b','FontSize', 12)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

%%%%%%%%%%%%%%%%%%%%%%%%


figure(9);
hold on
plot(x2(3,:),'b','LineWidth',2)
hold on
plot(V1(3,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12);
ylabel('Heave velocity [m.s^{-1}]','FontWeight','b','FontSize', 12);
legend('Actual', 'Virtual control','FontWeight','b','FontSize', 12)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.3 .26  .4 .3]);
plot(ax1,x2(3,:),'b','LineWidth',2);
grid off;box on
hold on
plot(ax1,V1(3,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 4000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);
hold on
plot(x2(4,:),'b', 'LineWidth',2)
hold on
plot(V1(4,:),'--r','LineWidth',2)
xlabel('Time [s]','FontWeight','b','FontSize', 12 );
ylabel('Pitch velocity [rad.s^{-1}]','FontWeight','b','FontSize', 12 );
legend('Actual', 'Virtual control','FontWeight','b','FontSize', 12 )
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.24 .39  .4 .3]);
plot(ax1,x2(4,:),'b','LineWidth',2);
grid off;box on
hold on
plot(ax1,V1(4,:),'--r','LineWidth',2);
ax1.XAxis.Exponent=0;
xlim([0 4000])
set( gca, 'FontWeight', 'b','FontSize', 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);
hold on
plot(x2(5,:),'b','LineWidth',2)
hold on
plot(V1(5,:),'--r','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('Yaw velocity [rad.s^{-1}]','FontWeight','b','FontSize', 12);
legend('Actual', 'Virtual control','FontWeight','b','FontSize', 12)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12);
hold on
subplot(3,1,1);
plot(u(1,:),'b','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_1 [N]','FontWeight','b','FontSize', 12);
legend('Control input','FontWeight','b','FontSize', 12)
xlim([0 NF])
ylim([-100 1000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

ax1 = axes('Position',[.2 .81  .4 .1]);
plot(ax1,u(1,:),'b','LineWidth',2);
grid off;box on
ax1.XAxis.Exponent=0;
xlim([5000 50000])
set( gca, 'FontWeight', 'b','FontSize', 12 );

subplot(3,1,2);
plot(u(2,:),'--r','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_2 [N.m]','FontWeight','b','FontSize', 12);
legend('Control input','FontWeight','b','FontSize', 12)
xlim([0 NF])
ylim([-1000 100])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

subplot(3,1,3);
plot(u(3,:),':g','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('\tau_3 [N.m]','FontWeight','b','FontSize', 12);
legend('Control input','FontWeight','b','FontSize', 12)
xlim([0 NF])
ylim([-100 1000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;
ax1 = axes('Position',[.2 .21  .4 .1]);
plot(ax1,u(3,:),':g','LineWidth',2);
grid off;box on
ax1.XAxis.Exponent=0;
xlim([5000 50000])
set( gca, 'FontWeight', 'b','FontSize', 12 );
%%%%%%%%%%%%%%%%%%%%%
figure(13);
hold on
plot(cost_com(:),'b','LineWidth',2)
xlabel('Time','FontWeight','b','FontSize', 12 );
ylabel('Cost','FontWeight','b','FontSize', 12);
legend('Cost','FontWeight','b','FontSize', 12)
xlim([0 NF])
set( gca, 'FontWeight', 'b','FontSize', 12 );
ax=gca;
ax.XAxis.Exponent=2;

% 
% for i=1:11
%     h=figure(i);
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperSize', [4 3]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 4 3]);
%     
%     set(h,'Units','Inches');
%     pos = get(h,'Position');
%     %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(h,"fig"+num2str(i),'-dpdf','-r300')
%      print(h,"fig"+num2str(i),'-dpng','-r300')
% end
% 
% h=figure(6);
% exportgraphics(h,'fig6.pdf' ,'Resolution',300)
% % %%%%%%%%%%%%%%%%%%
%     h=figure(12);
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperSize', [4 3]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 4 3]);
%     
%     set(h,'Units','Inches');
%     pos = get(h,'Position');
%     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(h,"fig12",'-dpdf','-r300')
%     print(h,"fig12",'-dpng','-r300')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure('WindowState','Maximized');
% grid on;
% xlabel('x', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
% ylabel('y', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
% zlabel('z', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
% hold on
% axis equal;
% 
% h_traj = plot3(eta(1,1), eta(2,1), eta(3,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% h_ref = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% view(57,32)
% % Plot the initial point of the AUV
% h_tr = animatedline('Color','b','LineStyle','--','LineWidth',3);
% h_re = animatedline('Color','r','LineStyle',':','LineWidth',3);
% set( gca, 'FontWeight', 'b','FontSize', 12 );
%  legend('AUV', 'Reference Point', 'Trajectory', 'Reference ', 'Location', 'northwest','FontWeight','b','FontSize',12 ,'Interpreter','tex');
% f = cell(NF/10,1) ; 
%  for t = 2:NF
%     if mod(t,100)==0
%     set(h_traj, 'XData', eta(1,t), 'YData', eta(2,t), 'ZData', eta(3,t));
%     addpoints(h_tr,eta(1,t),eta(2,t),eta(3,t));
% 
%     set(h_ref, 'XData', eta_d(1,t), 'YData', eta_d(2,t), 'ZData', eta_d(3,t));
%     addpoints(h_re,eta_d(1,t),eta_d(2,t),eta_d(3,t));
%     
%     drawnow
%     f{t/10} = getframe(gcf) ;
%     t
%     end
% end
% 
% obj = VideoWriter('AUV.avi');
% % obj.profile='Uncompressed AVI'
% obj.Quality = 100;
% obj.FrameRate = 50;
% open(obj);
% for i = 1:length(f)
%     writeVideo(obj, f{i}) ;
% end
% obj.close();
% 
% % figure(13);
% % grid on;
% % xlabel('x', 'Interpreter', 'latex');
% % ylabel('y', 'Interpreter', 'latex');
% % zlabel('z', 'Interpreter', 'latex');
% % axis equal;
% % hold on;
% % 
% % % Plot the full trajectory and reference path first
% % % plot3(eta(1,:), eta(2,:), eta(3,:), 'b:', 'LineWidth', 2);
% % % plot3(eta_d(1,:), eta_d(2,:), eta_d(3,:), 'r:', 'LineWidth', 2);
% % 
% % % Plot the initial point of the AUV
% % h_traj1 = plot3(eta(1,1), eta(2,1), eta(3,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% % h_ref1 = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% % 
% % 
% %  legend('Trajectory', 'Reference', 'AUV', 'Reference Point', 'Location', 'northwest');
% %  h_traj = plot3(eta(1,1), eta(2,1), eta(3,1), 'b');
% % h_ref = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1),  'r:');
% % 
% % % Move the points along the trajectory
% % for t = 2:NF
% %     % Update the position for the AUV and reference point
% % %     set(h_traj, 'XData', eta(1,1:t), 'YData', eta(2,1:t), 'ZData', eta(3,1:t));
% % %     set(h_ref, 'XData', eta_d(1,1:t), 'YData', eta_d(2,1:t), 'ZData', eta_d(3,1:t));
% % %   
% % %      set(h_traj1, 'XData', eta(1,t), 'YData', eta(2,t), 'ZData', eta(3,t));
% % %     set(h_ref1, 'XData', eta_d(1,t), 'YData', eta_d(2,t), 'ZData', eta_d(3,t));
% % 
% %     % Force MATLAB to render the plot and introduce a pause for the animation effect
% %     h_traj1 = plot3(eta(1,t), eta(2,t), eta(3,t), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% % h_ref1 = plot3(eta_d(1,t), eta_d(2,t), eta_d(3,t), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% % %     drawnow;
% % %     pause(.0000001);
% % %     t
% % end
% % 
% % 
% % figure(13);
% % grid on;
% % xlabel('x', 'Interpreter', 'latex');
% % ylabel('y', 'Interpreter', 'latex');
% % zlabel('z', 'Interpreter', 'latex');
% % axis equal;
% % hold on;
% % plot3(eta_d(1,:), eta_d(2,:), eta_d(3,:), 'r:', 'LineWidth', 2);
% % % Plot the initial point of the AUV
% % h_traj = plot3(eta(1,1), eta(2,1), eta(3,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% % h_ref = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% % 
% % for t = 2:NF
% %     % Update the position for the AUV and reference point
% %      set(h_traj, 'XData', eta(1,t), 'YData', eta(2,t), 'ZData', eta(3,t));
% % %     set(h_ref, 'XData', eta_d(1,t), 'YData', eta_d(2,t), 'ZData', eta_d(3,t));
% %    
% %     % Plot a single line segment for each point on the trajectory
% %     if mod(t,10)==0
% %         line('XData', eta(1,1:t), 'YData', eta(2,1:t), 'ZData', eta(3,1:t), 'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);
% %     end
% %         %line('XData', eta_d(1,t-1:t), 'YData', eta_d(2,t-1:t), 'ZData', eta_d(3,t-1:t), 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
% %    
% %     % Force MATLAB to render the plot and introduce a pause for the animation effect
% %     drawnow;
% % %     view(50,20)
% %     t
% % end
% % 
% % 
% % figure(13);
% % grid on;
% % xlabel('x', 'Interpreter', 'latex');
% % ylabel('y', 'Interpreter', 'latex');
% % zlabel('z', 'Interpreter', 'latex');
% % axis equal;
% % hold on;
% % plot3(eta_d(1,:), eta_d(2,:), eta_d(3,:), 'r:', 'LineWidth', 2);
% % % Plot the initial point of the AUV
% % h_traj = plot3(eta(1,1), eta(2,1), eta(3,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% % h_ref = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% % 
% % for t = 2:NF
% %     % Update the position for the AUV and reference point
% %      set(h_traj, 'XData', eta(1,t), 'YData', eta(2,t), 'ZData', eta(3,t));
% % %     set(h_ref, 'XData', eta_d(1,t), 'YData', eta_d(2,t), 'ZData', eta_d(3,t));
% %    
% %     % Plot a single line segment for each point on the trajectory
% %     if mod(t,10)==0
% %         line('XData', eta(1,1:t), 'YData', eta(2,1:t), 'ZData', eta(3,1:t), 'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);
% %     end
% %         %line('XData', eta_d(1,t-1:t), 'YData', eta_d(2,t-1:t), 'ZData', eta_d(3,t-1:t), 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
% %    
% %     % Force MATLAB to render the plot and introduce a pause for the animation effect
% %     drawnow;
% % %     view(50,20)
% %     t
% % end
% 
