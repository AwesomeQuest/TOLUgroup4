% skilaverkefni 4
% part 4
clc;clear all;
L = 5;
T = 30;
m = 100;
n = 100;
D = 0.01;
v = 0.1;
h = L / m;
k = T / n;

alpha = (k * D)/(h^2);
beta = (v * k)/(2*h);

A = diag(1 + 2 * alpha * ones(m-1, 1)) + ...
    diag((-alpha - beta) * ones(m-2, 1), -1) + ...
    diag((-alpha + beta) * ones(m-2, 1), 1);

% upphafsskilyrði
x = linspace(0, L, m+1);
w = zeros(m+1, n+1);
w(:, 1) = exp(-((x-1).^2)/D);

for j = 2:n+1
    b = w(2:m, j-1);
    w(2:m, j) = A \ b;
    w(1, j) = 0;
    w(m+1, j) = 0;
end

% Animation parameters
FPS = 30;                 % Frames per second for playback
n_frames = n + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames

% Generate frames
max_conc = max(w(:));
for j = 1:n_frames
    % Clear the figure
    clf;
    hold on;
    
    % Plot pollutant concentration
    plot(x, w(:, j), 'b-', 'LineWidth', 2);
    
    % Set axis limits and labels
    axis([0, L, 0, max_conc * 1.1]); % Adjust limits based on expected results
    xlabel('Position (m)');
    ylabel('Concentration (kg/m^3)');
    title(sprintf('Pollutant Distribution at t = %.2f min', (j-1) * k));
    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% Play animation
close all;
movie(gcf, F, 1, FPS);

% Save animation as video
v = VideoWriter('Pollutant_Animation1.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

[~, max_idx] = max(w(:, end));
disp(['Maximum concentration at t = 30 min: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);

%% part 5
clc;clear all;
L = 5;
T = 30;
m = 1000;
n = 1000;
D = 0.01;
v = 0.1;
h = L / m;
k = T / n;

alpha = (k * D)/(h^2);
beta = (v * k)/(2*h);

main_diag = (1 + 2 * alpha) * ones(m-1, 1);
lower_diag = (-alpha - beta) * ones(m-1, 1);
upper_diag = (-alpha + beta) * ones(m-1, 1);

A = spdiags([lower_diag main_diag upper_diag], [-1 0 1], m-1, m-1);

% upphafsskilyrði
x = linspace(0, L, m+1);
w = zeros(m+1, n+1);
w(:, 1) = exp(-((x-1).^2)/D);

for j = 2:n+1
    b = w(2:m, j-1);
    w(2:m, j) = A \ b;
    w(1, j) = 0;
    w(m+1, j) = 0;
end

% Animation parameters
FPS = 30;                 % Frames per second for playback
n_frames = n + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames

% Generate frames
max_conc = max(w(:));
for j = 1:n_frames
    % Clear the figure
    clf;
    hold on;
    
    % Plot pollutant concentration
    plot(x, w(:, j), 'b-', 'LineWidth', 2);
    
    % Set axis limits and labels
    axis([0, L, 0, max_conc * 1.1]); % Adjust limits based on expected results
    xlabel('Position (m)');
    ylabel('Concentration (kg/m^3)');
    title(sprintf('Pollutant Distribution at t = %.2f min', (j-1) * k));
    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% Play animation
close all;
movie(gcf, F, 1, FPS);

% Save animation as video
v = VideoWriter('Pollutant_Animation2.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

[~, max_idx] = max(w(:, end));
disp(['Maximum concentration at t = 30 min: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);

%% part 6
