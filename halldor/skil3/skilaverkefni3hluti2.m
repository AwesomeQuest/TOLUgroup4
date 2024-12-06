% skilaverkefni 3 hluti 2
% part 1
clc;clear all;


%% part 2
clc;clear all;
y0=[40;9];T=50;n=3200;

[tvec, yvec] = RK4solveht(y0,T,n);

plot(tvec, yvec)
title("Holling-Tanner")
xlabel("Tímibil")
ylabel("Fjöldi dýra")
legend("Prey", "Predators")

%% part 3
clc;clear all;
T=120;n=3200;

y0_1=[2;2];
y0_2=[15;10];
y0_3=[20;15];


[tvec1, yvec1] = RK4solveht(y0_1,T,n);
[tvec2, yvec2] = RK4solveht(y0_2,T,n);
[tvec3, yvec3] = RK4solveht(y0_3,T,n);

hold on
plot(yvec1(:,1), yvec1(:,2))
plot(yvec2(:,1), yvec2(:,2))
plot(yvec3(:,1), yvec3(:,2))
xlabel("Prey")
ylabel("Predators")
legend("Upphafsskilyrði [2;2]", "Upphafsskilyrði [15;10]", "Upphafsskilyrði [20;15]")

%% part 4
clc;clear all;
y0 = [40; 9];T = 50;n = 3200;

[tvec1, yvec1] = RK4solvelv(y0, T, n);
[tvec2, yvec2] = RK4solveht(y0, T, n);

subplot(2,1,1)
plot(tvec1, yvec1)
title("Lotka-Volterra")
legend("Prey","Predators")
subplot(2,1,2)
plot(tvec2, yvec2)
title("Holling-Tanner")
legend("Prey","Predators")