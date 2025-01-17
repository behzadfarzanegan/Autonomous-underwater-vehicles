
clear all;
clc;
close all
rng(1)

eta(1,1)=0.1;
eta(1,2)=eta(1,1);

eta(2,1)=0.1;
eta(2,2)=eta(2,1);

eta(3,1)=0.1;
eta(3,2)=eta(3,1);

eta(4,1)=0.1;
eta(4,2)=eta(4,1);

eta(5,1)=0.1;
eta(5,2)=eta(5,1);

eta(6,1)=1;
eta(6,2)=eta(6,1);

%%%%%%%%%%%%%%%%%%%%%%%%%
eta_1d(1)=0;
eta_2d(1)=0;
eta_3d(1)=0;
eta_4d(1)=0;
eta_5d(1)=0;
eta_6d(1)=1;

eta_1d(2)=0;
eta_2d(2)=0;
eta_3d(2)=0;
eta_4d(2)=0;
eta_5d(2)=0;
eta_6d(2)=1;
%%%%%%%%%%%%%%%%%%%%%
uu(1)=0;v(1)=0;r(1)=0;
uu(2)=0; v(2)=0; r(2)=0;
nue(:,1)=[0 0 0]';
nue(:,2)=[0 0 0]';

x(:,1) = [eta(1,1)-eta_1d(1),eta(2,1)-eta_2d(1),eta(3,1)-eta_3d(1),eta(4,1)-eta_4d(1),eta(5,1)-eta_5d(1),eta(6,1)-eta_6d(1)]'; %error
x(:,2) = [eta(1,1)-eta_1d(1),eta(2,1)-eta_2d(1),eta(3,1)-eta_3d(1),eta(4,1)-eta_4d(1),eta(5,1)-eta_5d(1),eta(6,1)-eta_6d(1)]'; 




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

Neuron_Num_a6=3;
Neuron_Num_c6=4;
Vc6=[rand(Neuron_Num_c6,2)];%no*1
Wc6(1,:)=rand(1,Neuron_Num_c6)*.01;
Wc6(2,:)=Wc6(1,:);
Wa6(1,:)=0.01*rand(1,Neuron_Num_a6);
Wa6(2,:)=Wa6(1,:);
r6=.1;q6=100;ro6=2.2;kappa6=10.7;a6=2;

T=0.01;
NF=5000;
gamma=0.9;
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

u6(1)= Wa6(1,:)*[x(6,1) x(6,1)^2 x(6,1)^3]'-1*x(6,1);
u6(2)= u6(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:NF
eta_1d(k)= sin(pi*k*T/5);
eta_1d(k+1)= sin(pi*(k+1)*T/5);

eta_2d(k)= sin(2*pi*k*T/5);
eta_2d(k+1)= sin(2*pi*(k+1)*T/5);

eta_3d(k)= 0.05*k*T+0.05;
eta_3d(k+1)= 0.05*(k+1)*T+0.05;

eta_4d(k)= 0;
eta_4d(k+1)= 0;

eta_5d(k)= 0;
eta_5d(k+1)= 0;

eta_6d(k)= atan(eta_2d(k)/eta_1d(k));
eta_6d(k+1)= atan(eta_2d(k+1)/eta_1d(k+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J1(k) = Wc1(k,:)*[ x(1,k)^2  x(1,k)^4 x(1,k)^6 x(1,k)^8]';

    del_J1=Wc1(k,:)*[2*x(1,k)  4*x(1,k)^3 6*x(1,k)^5 8*x(1,k)^7]';

    alpha1=-gamma/(2*r1)*del_J1';
  
    Yk_1 = q1*x(1,k-1)^2 +r1*u1(:,k-1)^2;
    delta_x1 = gamma*[x(1,k)^2  x(1,k)^4 x(1,k)^6 x(1,k)^8]' - [x(1,k-1)^2  x(1,k-1)^4 x(1,k-1)^6 x(1,k-1)^8]';
    EJK1 = Yk_1+Wc1(k,:)*delta_x1;
    
    Wc1(k+1,:) =Wc1(k,:)-ro1*(delta_x1'*EJK1)/(delta_x1'*delta_x1 +1);



     u1(k)= Wa1(k,:)*[x(1,k) x(1,k)^2 eta_1d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde1=  u1(k)-alpha1;
    sigma1a=[x(1,k) x(1,k)^2 x(1,k)^3];
    Sigma1a=sigma1a/(1+sigma1a*sigma1a');

     Wa1(k+1,:)= Wa1(k,:)-kappa1*Sigma1a*alpha_tilde1;

    x(1,k+1) = eta(1,k) +T*(u1(k))-eta_1d(k+1); %y_e
  
% 
    eta(1,k+1)=x(1,k+1)+eta_1d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J2(k) = Wc2(k,:)*[ x(2,k)^2  x(2,k)^4 x(2,k)^6 x(2,k)^8]';

    del_J2=Wc2(k,:)*[2*x(2,k)  4*x(2,k)^3 6*x(2,k)^5 8*x(2,k)^7]';

    alpha2=-gamma/(2*r2)*del_J2';
  
    Yk_2 = q2*x(2,k-1)^2 +r2*u2(:,k-1)^2;
    delta_x2 = gamma*[x(2,k)^2  x(2,k)^4 x(2,k)^6 x(2,k)^8]' - [x(2,k-1)^2  x(2,k-1)^4 x(2,k-1)^6 x(2,k-1)^8]';
    EJK2 = Yk_2+Wc2(k,:)*delta_x2;
    
    Wc2(k+1,:) =Wc2(k,:)-ro2*(delta_x2'*EJK2)/(delta_x2'*delta_x2 +1);



     u2(k)= Wa2(k,:)*[x(2,k) x(2,k)^2 eta_2d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde2=  u2(k)-alpha2;
    sigma2a=[x(2,k) x(2,k)^2 x(2,k)^3];
    Sigma2a=sigma2a/(1+sigma2a*sigma2a');

     Wa2(k+1,:)= Wa2(k,:)-kappa2*Sigma2a*alpha_tilde2;

    x(2,k+1) = eta(2,k) +T*(u2(k))-eta_2d(k+1); %y_e
  
% 
    eta(2,k+1)=x(2,k+1)+eta_2d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J3(k) = Wc3(k,:)*[ x(3,k)^2  x(3,k)^4 x(3,k)^6 x(3,k)^8]';

    del_J3=Wc3(k,:)*[2*x(3,k)  4*x(3,k)^3 6*x(3,k)^5 8*x(3,k)^7]';

    alpha3=-gamma/(2*r3)*del_J3';
  
    Yk_3 = q3*x(3,k-1)^2 +r3*u3(:,k-1)^2;
    delta_x3 = gamma*[x(3,k)^2  x(3,k)^4 x(3,k)^6 x(3,k)^8]' - [x(3,k-1)^2  x(3,k-1)^4 x(3,k-1)^6 x(3,k-1)^8]';
    EJK3 = Yk_3+Wc3(k,:)*delta_x3;
    
    Wc3(k+1,:) =Wc3(k,:)-ro3*(delta_x3'*EJK3)/(delta_x3'*delta_x3 +1);



     u3(k)= Wa3(k,:)*[x(3,k) x(3,k)^2 eta_3d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde3=  u3(k)-alpha3;
    sigma3a=[x(3,k) x(3,k)^2 x(3,k)^3];
    Sigma3a=sigma3a/(1+sigma3a*sigma3a');

     Wa3(k+1,:)= Wa3(k,:)-kappa3*Sigma3a*alpha_tilde3;

    x(3,k+1) = eta(3,k) +T*(u3(k))-eta_3d(k+1); %y_e
  
% 
    eta(3,k+1)=x(3,k+1)+eta_3d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J4(k) = Wc4(k,:)*[ x(4,k)^2  x(4,k)^4 x(4,k)^6 x(4,k)^8]';

    del_J4=Wc4(k,:)*[2*x(4,k)  4*x(4,k)^3 6*x(4,k)^5 8*x(4,k)^7]';

    alpha4=-gamma/(2*r4)*del_J4';
  
    Yk_4 = q4*x(4,k-1)^2 +r4*u4(:,k-1)^2;
    delta_x4 = gamma*[x(4,k)^2  x(4,k)^4 x(4,k)^6 x(4,k)^8]' - [x(4,k-1)^2  x(4,k-1)^4 x(4,k-1)^6 x(4,k-1)^8]';
    EJK4 = Yk_4+Wc4(k,:)*delta_x4;
    
    Wc4(k+1,:) =Wc4(k,:)-ro4*(delta_x4'*EJK4)/(delta_x4'*delta_x4 +1);



     u4(k)= Wa4(k,:)*[x(4,k) x(4,k)^2 eta_4d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde4=  u4(k)-alpha4;
    sigma4a=[x(4,k) x(4,k)^2 x(4,k)^3];
    Sigma4a=sigma4a/(1+sigma4a*sigma4a');

     Wa4(k+1,:)= Wa4(k,:)-kappa4*Sigma4a*alpha_tilde4;

    x(4,k+1) = eta(4,k) +T*(u4(k))-eta_4d(k+1); %y_e
  
% 
    eta(4,k+1)=x(4,k+1)+eta_4d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J5(k) = Wc5(k,:)*[ x(5,k)^2  x(5,k)^4 x(5,k)^6 x(5,k)^8]';

    del_J5=Wc5(k,:)*[2*x(5,k)  5*x(5,k)^3 6*x(5,k)^5 8*x(5,k)^7]';

    alpha5=-gamma/(2*r5)*del_J5';
  
    Yk_5 = q5*x(5,k-1)^2 +r5*u5(:,k-1)^2;
    delta_x5 = gamma*[x(5,k)^2  x(5,k)^4 x(5,k)^6 x(5,k)^8]' - [x(5,k-1)^2  x(5,k-1)^4 x(5,k-1)^6 x(5,k-1)^8]';
    EJK5 = Yk_5+Wc5(k,:)*delta_x5;
    
    Wc5(k+1,:) =Wc5(k,:)-ro5*(delta_x5'*EJK5)/(delta_x5'*delta_x5 +1);



     u5(k)= Wa5(k,:)*[x(5,k) x(5,k)^2 eta_5d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde5=  u5(k)-alpha5;
    sigma5a=[x(5,k) x(5,k)^2 x(5,k)^3];
    Sigma5a=sigma5a/(1+sigma5a*sigma5a');

     Wa5(k+1,:)= Wa5(k,:)-kappa5*Sigma5a*alpha_tilde5;

    x(5,k+1) = eta(5,k) +T*(u5(k))-eta_5d(k+1); %y_e
  
% 
    eta(5,k+1)=x(5,k+1)+eta_5d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J6(k) = Wc6(k,:)*[ x(6,k)^2  x(6,k)^4 x(6,k)^6 x(6,k)^8]';

    del_J6=Wc6(k,:)*[2*x(6,k)  6*x(6,k)^3 6*x(6,k)^5 8*x(6,k)^7]';

    alpha6=-gamma/(2*r6)*del_J6';
  
    Yk_6 = q6*x(6,k-1)^2 +r6*u6(:,k-1)^2;
    delta_x6 = gamma*[x(6,k)^2  x(6,k)^4 x(6,k)^6 x(6,k)^8]' - [x(6,k-1)^2  x(6,k-1)^4 x(6,k-1)^6 x(6,k-1)^8]';
    EJK6 = Yk_6+Wc6(k,:)*delta_x6;
    
    Wc6(k+1,:) =Wc6(k,:)-ro6*(delta_x6'*EJK6)/(delta_x6'*delta_x6 +1);



     u6(k)= Wa6(k,:)*[x(6,k) x(6,k)^2 eta_6d(k+1)]';%-eta(1,k)/T-1*x(1,k)+eta_1d(k+1)/T;
%       u1(k)=-5*x(1,k)-(sin(psi(1,k))*ud-cos(psi(1,k))*vd)-u3(k)*x(2,k);

    alpha_tilde6=  u6(k)-alpha6;
    sigma6a=[x(6,k) x(6,k)^2 x(6,k)^3];
    Sigma6a=sigma6a/(1+sigma6a*sigma6a');

     Wa6(k+1,:)= Wa6(k,:)-kappa6*Sigma6a*alpha_tilde6;

    x(6,k+1) = eta(6,k) +T*(u6(k))-eta_6d(k+1); %y_e
  
% 
    eta(6,k+1)=x(6,k+1)+eta_6d(k+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    R1=[cos(eta(6,k))*cos(eta(5,k)) -sin(eta(6,k))*cos(eta(4,k))+sin(eta(4,k))*sin(eta(5,k))*cos(eta(6,k)) sin(eta(6,k))*sin(eta(4,k))+sin(eta(5,k))*cos(eta(6,k))*cos(eta(4,k))
        sin(eta(6,k))*cos(eta(4,k)) cos(eta(6,k))*cos(eta(4,k))+sin(eta(4,k))*sin(eta(5,k))*sin(eta(6,k)) -cos(eta(6,k))*sin(eta(4,k))+sin(eta(5,k))*sin(eta(6,k))*cos(eta(4,k))
        -sin(eta(5,k))                   sin(eta(4,k))*cos(eta(5,k))                                          cos(eta(4,k)*cos(eta(4,k)))      ];

    V1(:,k)=R1'*[u1(k) u2(k) u3(k)]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

figure(1);hold on;
 plot(eta(1,1:NF),'LineWidth',2);
 hold on
 plot(eta_1d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);hold on;
 plot(eta(2,1:NF),'LineWidth',2);
 hold on
 plot(eta_2d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);hold on;
 plot(eta(3,1:NF),'LineWidth',2);
 hold on
 plot(eta_3d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(4);hold on;
 plot(eta(4,1:NF),'LineWidth',2);
 hold on
 plot(eta_4d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5);hold on;
 plot(eta(5,1:NF),'LineWidth',2);
 hold on
 plot(eta_5d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(6);hold on;
 plot(eta(6,1:NF),'LineWidth',2);
 hold on
 plot(eta_6d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);hold on;
 plot3(eta(1,1:NF),eta(2,1:NF),eta(3,1:NF),'LineWidth',2);
 hold on
 plot3(eta_1d(1:NF),eta_2d(1:NF),eta_3d(1:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(8);hold on;
 plot(V1(1,2:NF),'LineWidth',2);
 grid off; box on
set(gca,'GridLineStyle','--')
ax=gca;
 ax.XAxis.Exponent=3;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%