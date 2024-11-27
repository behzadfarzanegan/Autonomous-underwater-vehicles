function [output1,G] = NLinear_sys_NL_gamma_bah(x1,delta,r1,k,eta)
%----------------------Linear system Dynamics-------------------------------
T=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-Varying Exponential Stabilization of the Position and Attitude of an Underactuated Autonomous Underwater Vehicle
m11=215; m22=265; m33=265; m44=40; m55=80; m66=80; Delta=0.02;

d11=70; d22=100; d33= 100; d44= 30; d55=50; d66=50; 

W=1813;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=x1(1);v=x1(2);w=x1(3);q=x1(4);r=x1(5);

x=eta(1);y=eta(2);z=eta(3);
theta=eta(4);psi=eta(5);

X_prop=delta(1);
K_prop=delta(2);
dr=delta(3);



X_ext= m22/m11*v*r-m33/m11*w*q-d11/m11*u+0.05*sin(0.005*k);
Y_ext= -m11/m22*u*r - d22/m22*v;
Z_ext = m11/m33*u*q-d33/m33*w;
M_ext=(m33-m11)/m55*u*w-d55/m55*q-Delta/m55*sin(eta(4)*0);
N_ext=(m11-m22)/m66*u*v-d66/m66*r;

G1=[1/m11 0     0
    0     0     0
    0     0     0
    0     1/m55 0
    0     0     1/m66];

G = T.*G1;
F=T*[X_ext Y_ext Z_ext M_ext N_ext]';

xk1 = x1+F+G*delta-r1;
output1 = xk1; 

end
