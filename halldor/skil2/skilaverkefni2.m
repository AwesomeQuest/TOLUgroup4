% skilaverkefni 2
% part 1
clc;clear all;
T = (2:10) .*10 + 273.15;
mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;

A = [ones(size(mu)) T'.^(-1) T'.^(-2)];

beta = A\log(mu)

RMSE = norm(A*beta - log(mu))

norm(F(beta))

%% part 2
x0=zeros(3,1);tol=1e-5;
beta2 = gaussnewton(x0, tol)
norm(F(beta2))

%% part3
% least square method
T = (2:10) .*10 + 273.15;
mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;

B = [1./T' T' T'.^2];

beta3 = B\log(mu)

norm(F2(beta3))

RMSE2 = norm(B*beta2 - log(mu))

% gaussnewton method
betahat2 = gaussnewton2(x0, tol)
norm(F2(betahat2))

%% part 4
V0 = 0.01;
T0 = 20 + 273.15;
T1 = 80 + 273.15;
L = 0.1;
n = 2;
a = 0;
y = 0:0.01:L; 

K = V0./simpson(a,L,n);
v_data = zeros(size(y));

for i = 1:length(y)
    v_data(i) = K * simpson(a,y(i),n);
end

plot(y, v_data)
xlabel("Vegalengd (m)")
ylabel("Hraði vökva (m/s)")
title("V(y)")

%% part 5
V0 = 0.01;
T0 = 20 + 273.15;
T1 = 80 + 273.15;
L = 0.1;
n2 = [2,4,8,16,32,64,128,256];
n_real = 2^14;
a= 0;
y2 = L/2; 

v_real = V0 * (simpson(a,y2,n_real)/simpson(a,L,n_real));
%K = V0./simpson(a,L,n);
v_data2 = zeros(size(n2));

for i = 1:length(n2)
    v_data2(i) = V0.*(simpson(a,y2,n2(i))/simpson(a,L,n2(i)));
end

error = abs(v_data2 - v_real)

ratios = error(1:end-1)./error(2:end)

log2(ratios)

%% part 6
V0 = 0.01;tol = 0.5e-10;a = 0;L = 0.1;y2 = L/2;

nefnari = adapquad(a,y2,tol, "numerator");
teljari = adapquad(a,L,tol, "denominator");

V = V0 * (nefnari/teljari)
allintervals = adapquad()

%% part 7
x = linspace(0, 0.1,11); % x positions

time = linspace(0, 1, 100); % Time steps from 0 to 1 second

% Initialize the figure
figure;
hold on;
plot(v_data, x, '-o', 'DisplayName', 'Target Shape'); % Plot the target shape
h = plot(zeros(size(x)), x, '-r', 'LineWidth', 2, 'DisplayName', 'Bending Line'); % Vertical line at x=0
legend;
xlabel([' ']);
ylabel([' ']);
title(['T0 = 20 ' char(176) 'C and T1 = 80 ' char(176) 'C']); % Title
hold off;

% Animate the bending of the line
for t = 1:length(time)
    % Gradually interpolate the bending line to the target shape
    bending_shape = v_data * time(t); % Gradually progress from 0 (vertical) to x (curve)
    
    % Update the bending line
    set(h, 'XData', bending_shape, 'YData', x); % Adjust XData to reflect the bending shape, YData remains y
    
    pause(0.05); % Pause for animation effect
end

%% part 8
% part 4 gerður aftur
V0 = 0.01;
T0_2 = 80 + 273.15;
T1_2 = 20 + 273.15;
L = 0.1;
n = 2;
a = 0;
y = 0:0.01:L; 

K = V0./simpson2(a,L,n);
v_data3 = zeros(size(y));

for i = 1:length(y)
    v_data3(i) = K * simpson2(a,y(i),n); 
end

plot(y, v_data3)
xlabel("Bil (m)")
ylabel("Hraði vökva (m/s)")
title("V(y)")

% part 7 gerður aftur
x = linspace(0, 0.1,11); % x positions

time = linspace(0, 1, 100); % Time steps from 0 to 1 second

% Initialize the figure
figure;
hold on;
plot(v_data3, x, '-o', 'DisplayName', 'Target Shape'); % Plot the target shape
h = plot(zeros(size(x)), x, '-r', 'LineWidth', 2, 'DisplayName', 'Bending Line'); % Vertical line at x=0
legend;
xlabel([' ']);
ylabel([' ']);
title(['T0 = 80 ' char(176) 'C and T1 = 20 ' char(176) 'C']); % Title
hold off;

% Animate the bending of the line
for t = 1:length(time)
    % Gradually interpolate the bending line to the target shape
    bending_shape = v_data3 * time(t); % Gradually progress from 0 (vertical) to x (curve)
    
    % Update the bending line
    set(h, 'XData', bending_shape, 'YData', x); % Adjust XData to reflect the bending shape, YData remains y
    
    pause(0.05); % Pause for animation effect
end