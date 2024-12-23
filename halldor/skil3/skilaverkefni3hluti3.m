% skilaverkefni 3 hluti 3
% part 1
clc;clear all;
U0 = [0.8; 0.1; 8];
tspan = [0, 400];
T = 200;
n = 3200;

[t, U] = RK4solveholltann(U0, T, n);

X = U(:, 1);
Y = U(:, 2);
Z = U(:, 3);

figure;
plot(t, X, 'r', t, Y, 'g', t, Z, 'b');
xlabel('Time');
ylabel('Population');
legend('Plants (X)', 'Herbivores (Y)', 'Carnivores (Z)');
title('Time-dependent solutions of Holling-Tanner equations');

%% part 2
figure;
plot(X, Y, 'm');
xlabel('Plants (X)');
ylabel('Herbivores (Y)');
title('Phase plot: X vs Y');

figure;
plot(Y, Z, 'c');
xlabel('Herbivores (Y)');
ylabel('Carnivores (Z)');
title('Phase plot: Y vs Z');

% 3D Plot
figure;
plot3(X, Y, Z, 'k');
xlabel('Plants (X)');
ylabel('Herbivores (Y)');
zlabel('Carnivores (Z)');
title('3D Phase plot: X, Y, Z');

%% part 3
b1_values = linspace(0.5, 3, 5); % Range of b1 values
figure;
for i = 1:length(b1_values)
    % Update b1 in the function (modify holling_tanner or use an anonymous function)
    b1 = b1_values(i);
    [t, U] = ode45(@(t, U) holling_tanner_vary_b1(t, U, b1), tspan, U0);
    
    % Extract Z for plotting
    Z = U(:, 3);
    
    % Plot Z over time for each b1
    subplot(length(b1_values), 1, i);
    plot(t, Z);
    xlabel('Time');
    ylabel('Carnivores (Z)');
    title(['Z Dynamics for b1 = ', num2str(b1)]);
end
