% clc;close all;
N = 100;
M = 100;
T = 120;
L = 5;

[us1,ts1,xs1] = iterdiffsjalf(T,N,M);



% for i = 1:N+1
%     plot(xs1,us1(:,i))
%     text(2.5,0.5,num2str(ts1(i)))
%     ylim([0,1])
%     xlim([0,5])
%     pause(0.1)
% end


% plot(xs1,us1(:,round(N/6*6)))
[w,ts,xs] = iterdiffsjalf(T,N,M);
max_conc = max(w(:));
FPS = 30;                 % Frames per second for playback
n_frames = N + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames
x = linspace(0, L, M+1);
for j = 1:n_frames
    % Clear the figure
    clf;
    hold on;
    
    % Plot pollutant concentration
    plot(x, w(:, j), 'b-', 'LineWidth', 2);
    
    % Set axis limits and labels
    axis([0, L, 0, max_conc * 1.1]); % Adjust limits based on expected results
    xlabel('Staðsetning [m]');
    ylabel('Þéttleiki [kg/m^3]');
%     title(sprintf('Dreifing mengunarefnisins við t = %.2f mín', (j-1) * k));
    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% Play animation
% close all;
% movie(gcf, F, 1, FPS);


% Save animation as video
v = VideoWriter('Pollutant_Animation2.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

[~, max_idx] = max(w(:, end));
disp(['Maximum concentration at t = 30 s: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);
