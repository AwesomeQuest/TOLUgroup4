clc;close all;
T = 30;

N = 100;
M = 1000;
k = T/N;


[us1,ts,xs1] = iterdiffv1(T,N,M);
% [us2,ts,xs2] = iterdiffv2(T,N,M);
% [us3,ts,xs3] = iterdiffv1(T,10*N,M);

% plot(xs1,us1(:,round(N/6*6)))

% Animation parameters
FPS = 15;                 % Frames per second for playback
n_frames = N + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames

for i = 1:N+1
    clf;
    hold on
    plot(xs1,us1(:,i),LineWidth=2)
%     hold on
%     plot(xs2,us2(:,i))
%     plot(xs3,us3(:,10*i))
%     legend("original", "modified", "many iterations original")
    ylim([0,1])
%     pause(0.1)
    xlabel('Staðsetning [m]');
    ylabel('Þéttleiki [kg/m^3]');
    title(sprintf('Dreifing mengunarefnisins við t = %.2f s', (i-1) * k));
    grid on
    F(i) = getframe(hfig);
%     hold off
end
% close all;

% for i = 0:5
% %     figure;
%     j = -i + 5
%     plot(xs1,us1(:,round(N/6*j)+1),LineWidth=2)
%     ylim([0 1])
%     xlabel("Staðsetning [m]")
%     ylabel("Þéttleiki [kg/m^3]")
%     grid on
%     txt = "d4plot" + j + ".svg"
% %     saveas(gcf,txt)
%     
% end

% Play animation
close all;
movie(gcf, F, 1, FPS);

% Save animation as video
v = VideoWriter('Pollutant_Animation_d4.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);
