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
n = 10000; T=100;

% Upphafsgildi
y0_1 = [5.3;5.3];

[t,y] = RKsolver2(y0_1,n,T);

plot(t,y)

figure;
plot(y(1,:),y(2,:))

%% Daemi 3

% Breytum w = 0.3;
w = 1;

% Nyr stodugleikapunktur!
c2 = (sqrt(29) - 1)/2
y_eq3 = [c2;c2]; % y_eq1 = [7;0]; helst

[t,y] = RKsolver2(y0_1,n,T);

plot(t,y)

figure;
plot(y(1,:),y(2,:))

%% Daemi 4

% Thegar w = 0.3 er kerfid stodugt
