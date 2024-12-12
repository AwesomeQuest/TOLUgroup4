% clc;close all;
N = 200;
M = 200;
T = 100;

[us1,ts1,xs1] = iterdiffsjalf(T,N,M);
[us2,ts2,xs2] = iterdiffsjalf2(T,N,M);

%% Plot Singles

clf;
hold on;

% Plot pollutant concentration
plot(xs1,us1(:,round(N/6*6+1)), 'b-', 'LineWidth', 2)

% Set axis limits and labels
axis([0, 5, 0, 1]); % Adjust limits based on expected results
xlabel('Staðsetning [m]');
ylabel('Þéttleiki [kg/m^3]');
%     title(sprintf('Dreifing mengunarefnisins við t = %.2f mín', (j-1) * k));
grid on;

hold off;



%% Make Movie

L = 5;
[w,ts,xs] = iterdiffsjalf(T,N,M);
max_conc = 1;
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
v = VideoWriter('Pollutant_Animationsjalf1.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

[~, max_idx] = max(w(:, end));
disp(['Maximum concentration at t = 30 s: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);
