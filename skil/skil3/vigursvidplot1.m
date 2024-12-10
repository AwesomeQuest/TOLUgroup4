%% Hluti 1 - Daemi 1
clc; clear; close all;

% y1: prey
% y2: predator

% Fastar
% alpha: Vaxtarhradi bradarinnar ef enginn randyr eru til stadar.
% beta: Hlutfall af brad sem hverfur af voldum randyra.
% gamma: Vaxtarhradi randyrsins ef brað er til stadar.
% delta: Danartidni randyra.
alpha = 0.5; beta = 0.01; gamma = 0.005; delta = 0.2;
axis_max = 100;

% Buum til grid fyrir y1 og y2
skref = 5; % Breytilegt skref fyrir mynd
[y1, y2] = meshgrid(0:skref:axis_max, 0:skref:axis_max); 

% Finnum gildi a y1' og y2' fyrir hvert y1 og y2
dy1 = alpha.*y1 - beta.*y1.*y2;
dy2 = gamma.*y1.*y2 - delta.*y2;

% Normaliserum vigrana fyrir styttri vigra
magnitude = sqrt(dy1.^2 + dy2.^2);
dy1_norm = dy1 ./ magnitude;
dy2_norm = dy2 ./ magnitude;

% plotta vigursvið með lituðum örvum
figure;
hold on;

% Skala örvar fyrir skýrari mynd
arrow_scale = 3; % stærð örva
line_width = 2; % þykkt örva
num_colors = 256; % fjöldi lita
cmap = jet(num_colors); % Nota jet colormap (gott fyrir litaskiptingu)
colormap(cmap); % setja colormap fyrir current figure
magnitude_normalized = (magnitude - min(magnitude(:))) / (max(magnitude(:)) - min(magnitude(:))); 
color_indices = round(magnitude_normalized * (num_colors - 1)) + 1;

% Plotta örvar með tilsvarandi lit
for i = 1:size(y1, 1)
    for j = 1:size(y2, 2)
        % Teikna ör með lituðu skafti
        quiver(y1(i, j), y2(i, j), dy1_norm(i, j)*arrow_scale, dy2_norm(i, j)*arrow_scale, 0, ...
               'Color', cmap(color_indices(i, j), :), 'LineWidth', line_width, 'MaxHeadSize', 1.5);
    end
end

% Colorbar fyrir magnitude
colorbar;
caxis([min(magnitude(:)), max(magnitude(:))]); % Stillir colorbar mörk
xlabel('Bráð (y_1)');
ylabel('Rándýr (y_2)');
title('Vigursvið fyrir Lotka-Volterra líkan');
grid on;
axis([-5 axis_max+5 -5 axis_max+5]);
axis square;

% vista mynd
% saveas(gcf, 'LotkaVolterra_Vigursvid_Phase.png');