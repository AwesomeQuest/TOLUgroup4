clc;clear;close all;
T = 100;
N = 100;
M = N;
k = T/N;


[us1,ts1,xs1] = iterdiffd8(T,N,M);
% [us1,ts,xs1] = iterdiffd8org(T,N,M);
% plot(xs1,us1(:,1))

for i = 1:6
    figure;
    j = -i + 7
    plot(xs1,us1(:,round(N/6*j)+1),LineWidth=2)
    ylim([0 1])
    xlabel("Staðsetning [m]")
    ylabel("Þéttleiki [kg/m^3]")
    grid on
    txt = "d8plot" + j + ".svg"
    saveas(gcf,txt)
end


% plot3(xs1,ts1,us1);
%% LETS ANIMATE THIS MF

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

% Play animation
close all;
movie(gcf, F, 1, FPS);

% Save animation as video
v = VideoWriter('Pollutant_Animation_d8.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);