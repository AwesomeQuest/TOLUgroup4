%% Dæmi 1 - Lota 2
clc; clear; close all;

% Fastar
alpha = 0.5; % α: Vaxtarhraði báðarinnar ef enginn rándýr eru til staðar
beta = 0.01; % β: Hlutfall af bráð sem hverfur af völdum rándýra
gamma = 0.005; % γ: Vaxtarhraði rándýrsins ef bráð er til staðar
delta = 0.2; % δ: Dánartíðni rándýra
axis_max = 200;

% Grid fyrir y1: bráð og y2: rándýr
[y1, y2] = meshgrid(0:10:axis_max+30, 0:10:axis_max); % Grid range [0, 100] með skref 5

% Jöfnur - breytingar fyrir bráð og rándýr
dy1 = alpha*y1 - beta*y1.*y2; % (náttúrulegur vöxtur)-(fórnarlamb rándýrs)
dy2 = gamma*y1.*y2 - delta*y2; % (vöxtur vegna framboðs bráðar)-(náttúrulegur dauði)

% magnitute fyrir vektorana
magnitude = sqrt(dy1.^2 + dy2.^2);

% normalize vektorana fyrir örvar
dy1_norm = dy1 ./ magnitude;
dy2_norm = dy2 ./ magnitude;

% Veljum tima og skrefafjolda
T = 50;
n = 500;

% Upphafsgildi: 
y0_1 = [40; 9];

% Leysum kerfid med eulersolve falli
[tvec, yvec] = eulersolve(y0_1, n, T);

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
axis([-5 axis_max+35 -5 axis_max-10]);
axis square;

hold on

plot(yvec(1,:),yvec(2,:),LineWidth=1.5);

% vista mynd
% saveas(gcf, 'LotkaVolterra_Vigursvid_Phase1.svg');