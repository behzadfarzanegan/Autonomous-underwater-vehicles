figure('WindowState','Maximized');
grid on;
xlabel('x', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
ylabel('y', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
zlabel('z', 'FontWeight','b','FontSize', 12 ,'Interpreter','tex');
hold on
axis equal;

h_traj = plot3(eta(1,1), eta(2,1), eta(3,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
h_ref = plot3(eta_d(1,1), eta_d(2,1), eta_d(3,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
view(57,32)
% Plot the initial point of the AUV
h_tr = animatedline('Color','b','LineStyle','--','LineWidth',3);
h_re = animatedline('Color','r','LineStyle',':','LineWidth',3);
set( gca, 'FontWeight', 'b','FontSize', 12 );
 legend('AUV', 'Reference Point', 'Trajectory', 'Reference ', 'Location', 'northwest','FontWeight','b','FontSize',12 ,'Interpreter','tex');
f = cell(NF/10,1) ; 
 for t = 2:NF
    if mod(t,10)==0
    set(h_traj, 'XData', eta(1,t), 'YData', eta(2,t), 'ZData', eta(3,t));
    addpoints(h_tr,eta(1,t),eta(2,t),eta(3,t));

    set(h_ref, 'XData', eta_d(1,t), 'YData', eta_d(2,t), 'ZData', eta_d(3,t));
    addpoints(h_re,eta_d(1,t),eta_d(2,t),eta_d(3,t));
    
    drawnow
    f{t/10} = getframe(gcf) ;
    t
    end
end