%% Video demo of navigation solution
% Requirements: Run main prior
close all
h = figure;

v = VideoWriter('trajSim.avi');
VidoWriter.FrameRate = 60;
open(v)

% Initialization
plot(P_ref_INS(2,:),P_ref_INS(1,:),'r');hold on
plot(GPS_sim(2,1),GPS_sim(1,1),'xk')
plot(X_IN_GN(5,1),X_IN_GN(4,1),'b','LineWidth',1)
plot(P_INS(2,1),P_INS(1,1),'r')
axis([-700 700 -700 700])
legend('Reference','GPS','INS/GNSS','INS','Autoupdate','off')
title('Trajectory simulation','Interpreter','latex')
ylabel('North $(x_2)$','Interpreter','latex')
xlabel('East $(x_1)$','Interpreter','latex')

% Video settings
axis tight manual
ax = gca;
% ax.NextPlot = 'replaceChildren';

step = 20;
loops = N_between/step*(N_GNSS-1);
M(loops) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';

iter = 0;
iter_movie = 0;
prev_iter = [];
for i = 2:N_GNSS
    plot(GPS_sim(2,i-1),GPS_sim(1,i-1),'xk')
    axis([GPS_sim(2,i)-300 GPS_sim(2,i)+300 GPS_sim(1,i)-300 GPS_sim(1,i)+300])
    for j = 1:step:N_between
        iter = iter + step;
        iter_movie = iter_movie + 1;
        plot(P_INS(2,iter),P_INS(1,iter),'.r','LineWidth',1)
        plot(X_IN_GN(5,iter),X_IN_GN(4,iter),'.b','LineWidth',1)

%         drawnow
%         M(iter_movie) = getframe;
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

end
axis([-700 700 -700 700])
for k = 1:60
    frame = getframe(gcf);
    writeVideo(v,frame);
end

% h.Visible = 'on'
% movie(M)

close(v)
beep
disp('Done')
% plot(P_ref_INS(2,:),P_ref_INS(1,:),'r');hold on
% plot(GPS_sim(2,:),GPS_sim(1,:),'.k')
% plot(X_IN_GN(5,:),X_IN_GN(4,:),'g','LineWidth',1)
% plot(P_INS(2,:),P_INS(1,:))
% plot(X_IN_GN_iter(5,:),X_IN_GN_iter(4,:),'.b','LineWidth',1)
% axis([-700 700 -700 700])
% legend('Reference','GPS','INS/GNSS','INS')
% title('Trajectory simulation')