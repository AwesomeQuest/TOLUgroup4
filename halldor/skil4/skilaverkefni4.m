% skilaverkefni 4
% part 4
clc;clear all;
L = 5;T = 30;m = 100;n = 100;
D = 0.01;v = 0.1;h = L / m;k = T / n;

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
    title(sprintf('Pollutant Distribution at t = %.2f s', (j-1) * k));
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
disp(['Maximum concentration at t = 30 s: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);

%% part 5
clc;clear all;
L = 5;T = 30;m = 1000;n = 1000;
D = 0.01;v = 0.1;h = L / m;k = T / n;

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
    title(sprintf('Pollutant Distribution at t = %.2f s', (j-1) * k));
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
disp(['Maximum concentration at t = 30 s: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);

%%
% Common parameters
clc; clear all;
L = 5; T = 30; 
D = 0.01; v = 0.1; 
time_steps = [0, T/6, 2*T/6, 3*T/6, 4*T/6, 5*T/6]; % Time points for plotting

% Parameters for Part 4
m1 = 100; n1 = 100;
h1 = L / m1; k1 = T / n1;
alpha1 = (k1 * D) / (h1^2);
beta1 = (v * k1) / (2 * h1);

A1 = diag(1 + 2 * alpha1 * ones(m1-1, 1)) + ...
    diag((-alpha1 - beta1) * ones(m1-2, 1), -1) + ...
    diag((-alpha1 + beta1) * ones(m1-2, 1), 1);

x1 = linspace(0, L, m1+1);
w1 = zeros(m1+1, n1+1);
w1(:, 1) = exp(-((x1-1).^2)/D);

for j = 2:n1+1
    b1 = w1(2:m1, j-1);
    w1(2:m1, j) = A1 \ b1;
    w1(1, j) = 0;
    w1(m1+1, j) = 0;
end

% Parameters for Part 5
m2 = 1000; n2 = 1000;
h2 = L / m2; k2 = T / n2;
alpha2 = (k2 * D) / (h2^2);
beta2 = (v * k2) / (2 * h2);

main_diag = (1 + 2 * alpha2) * ones(m2-1, 1);
lower_diag = (-alpha2 - beta2) * ones(m2-1, 1);
upper_diag = (-alpha2 + beta2) * ones(m2-1, 1);

A2 = spdiags([lower_diag main_diag upper_diag], [-1 0 1], m2-1, m2-1);

x2 = linspace(0, L, m2+1);
w2 = zeros(m2+1, n2+1);
w2(:, 1) = exp(-((x2-1).^2)/D);

for j = 2:n2+1
    b2 = w2(2:m2, j-1);
    w2(2:m2, j) = A2 \ b2;
    w2(1, j) = 0;
    w2(m2+1, j) = 0;
end

% Plot results for comparison
for t_idx = 1:length(time_steps)
    t_plot = time_steps(t_idx);
    idx1 = round(t_plot / k1) + 1; % Find index for Part 4
    idx2 = round(t_plot / k2) + 1; % Find index for Part 5
    
    figure;
    hold on;
    plot(x1, w1(:, idx1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dæmi 4 (m=n=100)');
    plot(x2, w2(:, idx2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Dæmi 5 (m=n=1000)');
    hold off;
    
    xlabel('Staðsetning (m)');
    ylabel('Þéttleiki (kg/m^3)');
    title(sprintf('Dreifing mengunarefninsin við t = %.0f/6T', t_idx-1));
    legend('Location', 'best');
    grid on;
end
