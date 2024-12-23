%% Hluti 2 - Daemi 1
clc; clear; close all;

% y1: prey
% y2: predator

% Fastar
r_1 = 1; r_2 = 0.1; k = 7;
d = 1; j = 1; w = 0.3;

% Fundum stodugleikapunkta a bladi
c1 = -(-39-sqrt(4321))/20;
y_eq1 = [k;0]; y_eq2 = [c1;c1];

%% Daemi 2

% Veljum skrefafjolda og lokatima
n = 600; T=120;

% Upphafsgildi
y0_1 = [10; 10]; y0_2 = [4; 4]; y0_3 = [2; 2];

[t,y1] = RKsolver2(y0_1,n,T);
[t,y2] = RKsolver2(y0_2,n,T);
[t,y3] = RKsolver2(y0_3,n,T);

figure;
hold on
plot(t,y1(1,:),'b', t, y1(2,:), '--b', 'LineWidth', 1.5);
plot(t,y2(1,:),'g', t, y2(2,:), '--g', 'LineWidth', 1.5);
plot(t,y3(1,:),'r', t, y3(2,:), '--r', 'LineWidth', 1.5);
xlim([0 45]);

legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2');
xlabel('Tími');
ylabel('Stofn');
title("Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")

% figure;
% plot(y1(1,:),y1(2,:))

%% Daemi 3

% Breytum w = 0.3;
w = 1;

% Nyr stodugleikapunktur!
c2 = (sqrt(29) - 1)/2
y_eq3 = [c2;c2]; % y_eq1 = [7;0]; helst
[t,y1] = RKsolver2(y0_1,n,T,w);
[t,y2] = RKsolver2(y0_2,n,T,w);
[t,y3] = RKsolver2(y0_3,n,T,w);

figure;
hold on
plot(t,y1(1,:),'b', t, y1(2,:), '--b', 'LineWidth', 1.5);
plot(t,y2(1,:),'g', t, y2(2,:), '--g', 'LineWidth', 1.5);
plot(t,y3(1,:),'r', t, y3(2,:), '--r', 'LineWidth', 1.5);

legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2');
xlabel('Tími');
ylabel('Stofn');
title("Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")

figure;
hold on
plot(y1(1,:),y1(2,:),'b',LineWidth=1.1)
plot(y2(1,:),y2(2,:),'g',LineWidth=1.1)
plot(y3(1,:),y3(2,:),'r',LineWidth=1.1)

legend("y_0=[10 10]", "y_0=[4 4]", "y_0=[2 2]")
xlabel("Bráð");
ylabel('Rándýr');
title("Phase plot Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")


%% Daemi 4

n=3200;T=200;
[t, rk_1]=RKsolver(y0_1,n,T);
[time2, rk_2]=RKsolver(y0_2,n,T);
[time3, rk_3]=RKsolver(y0_3,n,T);

figure;
hold on
plot(t, rk_1(1, :), 'b', t, rk_1(2, :), '--b', 'LineWidth', 1.5);
plot(t, rk_2(1, :), 'g', t, rk_2(2, :), '--g', 'LineWidth', 1.5); % dæmi 2
plot(t, rk_3(1, :), 'r', t, rk_3(2, :), '--r', 'LineWidth', 1.5); % dæmi 3
legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2','location','north');
xlabel('Tími');
ylabel('Stofn');
title("Lotka-Volterra líkan leyst með aðferð Runge-Kutta (RK4)")
xlim([0 60])
hold off
figure;
plot(rk_1(1,:),rk_1(2,:),'b',LineWidth=1.1)
hold on
plot(rk_2(1,:),rk_2(2,:),'g',LineWidth=1.1)
plot(rk_3(1,:),rk_3(2,:),'r',LineWidth=1.1)
legend("y_0=[10 10]", "y_0=[4 4]", "y_0=[2 2]")
xlabel("Bráð");
ylabel('Rándýr');
title("Phase plot Lotka-Volterra líkan leyst með aðferð Runge-Kutta (RK4)")
hold off
