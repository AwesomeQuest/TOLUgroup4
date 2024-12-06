% skilaverkefni 3 hluti 1
% part 1
clc;clear all;
alpha=0.5;beta=0.01;gamma=0.005;delta=0.2;
[y1, y2] = meshgrid(0:10:300, 0:10:300);

dy1 = alpha.*y1 - beta.*y1.*y2;
dy2 = gamma.*y1.*y2 - delta.*y2;

% normalizeum vigrana
lengd_dy = sqrt(dy1.^2 + dy2.^2);
dy1_norm = dy1./lengd_dy;
dy2_norm = dy2./lengd_dy;

figure;
quiver(y1, y2, dy1_norm, dy2_norm)
xlabel("Prey (y1)")
ylabel("Predators (y2)")

%% part 2
y0 = [40; 9];T = 50;n = 500;

[tvec, yvec] = eulersolve(y0, n, T);

plot(tvec, yvec)
title("Lotka-Volterra")
xlabel("Tímibil")
ylabel("Fjöldi dýra")
legend("Prey", "Predators")

%% part 3
y0 = [40; 9];T = 50;n = 500;

[tvec2, yvec2] = eulersolve(y0, n, T);

hold on
plot(yvec2(1, :), yvec2(2, :));
quiver(y1, y2, dy1_norm, dy2_norm);
xlabel("Prey")
ylabel("Predators")

%% part4
clc;clear all;
y0 = [40; 9]; 
T = 50; 
%n = [100, 200, 400, 800];
n = [6400, 12800, 25600, 51200];
n_real = 10^5; 

[tref, yref] = RK4solvelv(y0, T, n_real);  

error = zeros(1, length(n));

for i = 1:length(n)
    [tvec3, yvec3] = eulersolve(y0, n(i), T);

    yref_interp = interp1(tref, yref, tvec3);
    yref_interp = yref_interp';

    %error(i) = max(abs(yvec3(:,end) - yref_interp(:,end)));
    error(i) = norm(yvec3(:,end) - yref_interp(:,end));
end

disp('Errors:');
disp(error);

ratio = error(1:end-1) ./ error(2:end);
disp('Error ratios:');
disp(ratio);

% plottum up errorið
h_values = T ./ n; 
figure;
loglog(h_values, error, '-o');
xlabel('Step size (h)');
ylabel('Error');
title('Error Convergence of Euler Method');

% til þess að sjá hvaða stig þessi aðferð er
p = polyfit(log(h_values), log(error), 1); 
fprintf('Convergence: %.2f\n', p(1));

%% part 5
clc;clear all;
y0 = [40; 9];T = 50;n = 500;n2=10^4;

y = RKsolver(y0,T,n);
t_rk4 = linspace(0,T,n+1);
[tvec, yvec] = eulersolve(y0,n2,T);

hold on
plot(t_rk4, y,"--")
plot(tvec, yvec)
legend("Prey RK4", "Predators RK4", "Prey euler", "Predators euler")

%% part 6
clc;clear all;
y0 = [40; 9]; 
T = 50; 
n = [100, 200, 400, 800]; 
n_ref = 10^6; 
tspan = [0,T];

[tref, yref] = RK4solvelv2(y0,n_ref,T);
%[tref, yref] = ode45(@lv,tspan,y0);

error = zeros(1, length(n));

for i = 1:length(n)
    [tvec, yvec] = RK4solvelv2(y0,n(i),T);

    y_ref_interp = interp1(tref, yref', tvec)';
    
    %error(i) = norm(yvec(:,end) - y_ref_interp(:,end));
    error(i) = max(abs(yvec(:,end) - y_ref_interp(:,end))); 
end


disp('Errors:');
disp(error);

ratios = error(1:end-1) ./ error(2:end);
disp('Error ratios:');
disp(ratios);

% plottum up errorið
h_values = T ./ n;
figure;
loglog(h_values, error, '-o');
xlabel('Step size (h)');
ylabel('Error');
title('Error Convergence for RK4');

% til þess að sjá hvaða stig þessi aðferð er
p = polyfit(log(h_values), log(error), 1);
fprintf('Convergence rate: %.2f\n', p(1));


