function [eta] = kinematic(eta,u)
T=0.01;
   R1=[cos(eta(5))*cos(eta(4)) -sin(eta(5)) sin(eta(4))*cos(eta(5))
        sin(eta(5))*cos(eta(4)) cos(eta(5)) sin(eta(4))*sin(eta(5))
        -sin(eta(4))              0             cos(eta(4))];

    R2=[1 0
        0 1/cos(eta(4))];

        %%%%%%%%%%%%%%%%%

     eta(1:3)=eta(1:3)+T*R1*u(1:3);
     eta(4:5)=eta(4:5)+T*R2*u(4:5);
    
