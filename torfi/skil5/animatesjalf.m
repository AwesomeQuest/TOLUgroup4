% clc;close all;


%% Plot Singles v1
N = 200;
M = 200;
T = 100;
L = 5;
[us1,ts1,xs1] = iterdiffsjalf(T,N,M);

clf;
hold on;

% Plot pollutant concentration
plot(xs1,us1(:,round(N/6*6+1)), 'b-', 'LineWidth', 2)

% Set axis limits and labels
axis([0, 5, 0, 1]); % Adjust limits based on expected results
xlabel('Staðsetning [m]');
ylabel('Þéttleiki [kg/m^3]');

grid on;

hold off;


%% Plot Singles v2
N = 200;
M = 200;
T = 100;
L = 5;

[us2,ts2,xs2] = iterdiffsjalf2(T,N,M);
clf;
hold on;

% Plot pollutant concentration
plot(xs2,us2(:,round(N/6*6+1)), 'b-', 'LineWidth', 2)

% Set axis limits and labels
axis([0, 5, 0, 1]); % Adjust limits based on expected results
xlabel('Staðsetning [m]');
ylabel('Þéttleiki [kg/m^3]');

grid on;

hold off;



%% Make Movie v1
N = 200;
M = 200;
T = 100;
L = 5;

[w,ts,xs] = iterdiffsjalf(T,N,M);
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
    axis([0, L, 0, 1]); % Adjust limits based on expected results
    xlabel('Staðsetning [m]');
    ylabel('Þéttleiki [kg/m^3]');

    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% % Play animation
% close all;
% movie(gcf, F, 1, FPS);


% Save animation as video
v = VideoWriter('Pollutant_Animationsjalf1.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

%% Make Movie v2
N = 200;
M = 200;
T = 100;
L = 5;

[w,ts,xs] = iterdiffsjalf2(T,N,M);
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
    axis([0, L, 0, 1]); % Adjust limits based on expected results
    xlabel('Staðsetning [m]');
    ylabel('Þéttleiki [kg/m^3]');

    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% % Play animation
% close all;
% movie(gcf, F, 1, FPS);


% Save animation as video
v = VideoWriter('Pollutant_Animationsjalf2.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);
