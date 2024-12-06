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
axis_max = 200;

% Buum til grid fyrir y1 og y2
skref = 5; % Breytilegt skref fyrir mynd
[y1, y2] = meshgrid(0:skref:axis_max, 0:skref:axis_max); 

% Finnum gildi a y1' og y2' fyrir hvert y1 og y2
y1_diff = alpha.*y1 - beta.*y1.*y2;
y2_diff = gamma.*y1.*y2 - delta.*y2;

% Normaliserum vigrana fyrir styttri vigra
norm_factor = sqrt(y1_diff.^2 + y2_diff.^2);
y1_diff_normalized = y1_diff ./ norm_factor;
y2_diff_normalized = y2_diff ./ norm_factor;

% Plottum vigursvidid
figure;
quiver(y1, y2, y1_diff_normalized, y2_diff_normalized, 'b');
% quiver(y1, y2, y1_diff, y2_diff, 'b');
xlabel('y_1');
ylabel('y_2');
title('Vigursvið diffurjafnanna');
axis([0 axis_max 0 axis_max]);
grid on;

%% Daemi 2

% Veljum tima og skrefafjolda
T = 50;
n = 500;

% Upphafsgildi: 
y0_1 = [40; 9];

% Leysum kerfid med eulersolve falli
[tvec1, yvec1] = eulersolve(y0_1, n, T);

% Plottum nidurstodur
figure;
plot(tvec1,yvec1(1,:),'b',LineWidth=1.5);
hold on
plot(tvec1,yvec1(2,:),'r',LineWidth=1.5);
% legend("y_1 (með y_0_1 = 40)","y_2 (með y_0_2 = 9)","y_1 (með y_0_1 = 9)","y_2 (með y_0_2 = 40)")
xlabel('Tími');
ylabel('Fjöldi');
% title('Lausn á upphafsgildisverkefni');
xlim([0 T]);
grid on

% Profum nytt upphafsgildi
y0_2 = [25; 25];
[tvec2, yvec2] = eulersolve(y0_2, n, T);
hold on
plot(tvec2,yvec2(1,:),'--g',LineWidth=1.5);
plot(tvec2,yvec2(2,:),'--m',LineWidth=1.5);
legend("y_1 (með y_0_1 = 40)","y_2 (með y_0_2 = 9)","y_1 (með y_0_1 = 25)","y_2 (með y_0_2 = 25)")

% saveas(gcf,"hluti1_d2_plot.png");

figure;
plot(yvec1(1,:),yvec1(2,:),'b',LineWidth=1.5)
hold on
plot(yvec2(1,:),yvec2(2,:),'r',LineWidth=1.5)
xlabel("y_1")
ylabel('y_2')
title('Fasaferilsmynd með aðferð Eulers')
legend("y_0 = [40; 9]","y_0 = [25; 25]")
grid on

saveas(gcf,"hluti1_d2_plot2.png");

% Plottum nidurstodur
% figure;
% plot(tvec2,yvec2,LineWidth=1.5);
% legend("y_1 (Bráð)","y_2 (Rándýr)")
% xlabel('Tími');
% ylabel('Fjöldi');
% title('Lausn á upphafsgildisverkefni');
% xlim([0 T]);
% grid on

%% Daemi 3

figure;
quiver(y1, y2, y1_diff_normalized, y2_diff_normalized, 'b');
hold on
plot(yvec1(1,:),yvec1(2,:),LineWidth=1.5);
% hold on
% plot(y_ref(:,1),y_ref(:,2),LineWidth=1.5);

xlabel('y_1');
ylabel('y_2');
axis([0 axis_max 0 axis_max]);

%% Daemi 4

% Faum viðmiðunargildi til thess ad bera saman (finna skekkju)
tspan = [0 T];
[t_ref, y_ref] = RKsolver(y0_1,100,T);

% Buum til vigur med mismunandi gildum a n
% n_vector = [100 200 400 800];
n_vector = 100*2.^(0:12);

% Tomt fylki til ad geyma mismuninn
% error1 = zeros(size(n_vector));
error1 = zeros(length(y0_1),length(n_vector));
% index1 = find(abs(t_ref - 25) < 1e-10);

for i = 1:length(n_vector)
    [t_vec, y_vec] = eulersolve(y0_1, n_vector(i), T);
%     index = find(t_vec == 25);
%     error_vec(:,i) = sum(abs(y_ref(:,end) - y_vec(:,end)));
%     error1(i) = norm(y_ref(:,end) - y_vec(:,end));
    error1(:,i) = abs(y_ref(:,end) - y_vec(:,end));
end

% index = find(t_vec == 30)

factor_munur = error1(:,1:end-1) ./ error1(:,2:end);
log2_ratio = log2(factor_munur);

% plot(t_vec,y_vec);
% hold on
% plot(t_ref,y_ref);
% 
% legend("y_1 Euler","y_2 Euler","y_1 RK4","y_2 RK4");

% plot(n_vector,error_vec);

%% Daemi 5

T1 = 100;
y0_loop = [y0_1 y0_2 [100; 25]];

for i = 1:length(y0_loop)
    [t_eul, y_eul] = eulersolve(y0_loop(:,i),16000,T1);
    [t_RK, y_RK] = RKsolver(y0_loop(:,i),1000,T1);
    figure;
    plot(t_eul,y_eul(1,:),'--b',LineWidth=1.5)
    hold on
    plot(t_eul,y_eul(2,:),'--r',LineWidth=1.5)
    plot(t_RK,y_RK(1,:),'g',LineWidth=1.5)
    plot(t_RK,y_RK(2,:),'m',LineWidth=1.5)
    legend("y_1 Euler","y_2 Euler","y_1 RK4","y_2 RK4");
    title("Upphafsgildi: y_1 = " + y0_loop(1,i) + " y_2 = " + y0_loop(2,i))
    xlabel("Tími")
    ylabel("Fjöldi")
    xlim([0 T1+1])
%     figure;
%     plot(y_eul(1,:),y_eul(2,:),'',LineWidth=1.5)
%     hold on
%     plot(y_RK(1,:),y_RK(2,:),'',LineWidth=1.5)
%     title("Upphafsgildi: y_1 = " + y0_loop(1,i) + " y_2 = " + y0_loop(2,i))
end

% saveas(gcf,'hluti1_d5_plot.png')

%% Daemi 6

% Tomt fylki til ad geyma mismuninn
% error_RK = zeros(size(n_vector));
error_RK = zeros(length(y0_1),length(n_vector));

for i = 1:length(n_vector)
    [t_vec, y_vec] = RKsolver(y0_1, n_vector(i), T);
%     error_vec(:,i) = sum(abs(y_ref(:,end) - y_vec(:,end)));
%     disp(y_vec(:,end));
%     error_RK(i) = norm(y_ref(:,end) - y_vec(:,end));
    error_RK(:,i) = abs(y_ref(:,end) - y_vec(:,end));
end

factor_RK = error_RK(:,1:end-1) ./ error_RK(:,2:end);
log2_ratio_RK = log2(factor_RK);

% 4. stigs baeting thar til 100*2^7