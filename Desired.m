function [output1,output2] = Desired(eta, x1,k)
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


%%%%%%%%%%%%%%%%%%%%
xd=x+T*(cos(psi)*cos(theta)*u-sin(psi)*v+sin(theta)*cos(psi)*w);
yd=y+T*(sin(psi)*cos(theta)*u+cos(psi)*v+sin(theta)*sin(psi)*w);
zd=z+T*(-sin(theta)*u+cos(theta)*w);

thetad=theta+T*(q);
psid=psi+T*(r/cos(theta));

output1=[xd yd zd thetad psid]';

%%%%%%%%%%%%%%%%%%


X_ext= m22/m11*v*r-m33/m11*w*q-d11/m11*u;
Y_ext= -m11/m22*u*r - d22/m22*v;
Z_ext = m11/m33*u*q-d33/m33*w;
M_ext=(m33-m11)/m55*u*w-d55/m55*q-Delta/m55*sin(eta(4));
N_ext=(m11-m22)/m66*u*v-d66/m66*r;

G1=[1/m11 0     0
    0     0     0
    0     0     0
    0     1/m55 0
    0     0     1/m66];

G = T.*G1;
F=T*[X_ext Y_ext Z_ext M_ext N_ext]';

if k<=20000

tau=[.5*d11-(m22*v*r-m33*w*q) -(m33-m11)*u*w-d55*q-69.42*sin(theta)+m55*(-theta +.2-(thetad-theta)/T) -(m11-m22)*u*v+15*sin(0.061*T*k)-15*cos(0.061*k*T)]';

elseif k> 20000 && k<=30000

tau=[.5*d11-(m22*v*r-m33*w*q) 0 0]';
else 
tau=[.5*d11-(m22*v*r-m33*w*q) -(m33-m11)*u*w-d55*q-69.42*sin(theta)+m55*(-theta +.2-(thetad-theta)/T) -(m11-m22)*u*v+15*sin(0.061*T*k)-15*cos(0.061*k*T)]';

end
output2 = x1+F+G*tau;


end
