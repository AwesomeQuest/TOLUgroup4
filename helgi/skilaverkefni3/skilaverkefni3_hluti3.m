%% Hluti 3 - Daemi 1
clc; clear; close all

% X: plontur
% Y: plontuaetur
% Z: kjotaetur

% Fastar:
a_1 = 5; b_1 = 3; a_2 = 0.1;
b_2 = 2; d_1 = 0.4; d_2 = 0.01;

% Gefum okkur lokatima, skrefastaerd og upphafsgildi
T = 400; n = 100000; y0_1 = [0.8; 0.1; 8];

% Leysum med RK4
[t1,y1] = RKsolver3(y0_1,n,T);

% Plottum graf breytanna sem fall af tima
plot(t1,y1, LineWidth=1.5)
xlabel("Tími [ótilgreind eining]")
ylabel("Fjöldi")
xlim([0 T+1])
title("Holling-Tanner líkanið - RK4")
legend("Plöntur (X)", "Plöntuætur (Y)", "Kjötætur (Z)")

saveas(gcf,"h3_d1_plot.png")

%% Daemi 2

% Skilgreinum ny upphafsgildi
y0_2 = [2; 0.5; 12];
y0_3 = [0.1; 0.1; 5];

% Leysum med RK4
[t2, y2] = RKsolver3(y0_2,n,T);
[t3, y3] = RKsolver3(y0_3,n,T);

% Teiknum fasaferla 3ja lausna
figure;
plot3(y1(1,:),y1(2,:),y1(3,:),LineWidth=1.5)
xlabel("X")
ylabel("Y")
zlabel("Z")
% title("Upphafsgildi: y_1 = [" + y0_1(1) + "; " + y0_1(2) + "; " + y0_1(3) + "]")
saveas(gcf,"d2_fasaferill1.png")
figure;
plot3(y2(1,:),y2(2,:),y2(3,:),LineWidth=1.5)
xlabel("X")
ylabel("Y")
zlabel("Z")
% title("Upphafsgildi: y_2 = [" + y0_2(1) + "; " + y0_2(2) + "; " + y0_2(3) + "]")
saveas(gcf,"d2_fasaferill2.png")
figure;
plot3(y3(1,:),y3(2,:),y3(3,:),LineWidth=1.5)
xlabel("X")
ylabel("Y")
zlabel("Z")
% title("Upphafsgildi: y_3 = [" + y0_3(1) + "; " + y0_3(2) + "; " + y0_3(3) + "]")
saveas(gcf,"d2_fasaferill3.png")

%% Daemi 3

b_1_vec = 0.5:0.5:3;

for i = 1:length(b_1_vec)
    [t, y] = RKsolver3(y0_1,n,T,b_1_vec(i));
    figure;
    plot(t,y,LineWidth=1.5)
    xlim([0 T+1])
    title("b_1 = " + b_1_vec(i))
    figure;
    plot3(y(1,:),y(2,:),y(3,:))
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    title("b_1 = " + b_1_vec(i))
end
