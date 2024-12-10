clc; clear; close all;

[t,y1] = RK4solveCLV([0.1;0.2;0.2], 1000, 200.0);

plot(t,y1(1,:),'b', t, y1(2,:), 'r',t,y1(3,:),'g', 'LineWidth', 1.5)
legend('Tegund 1', 'Tegund 2', 'Tegund 3');

