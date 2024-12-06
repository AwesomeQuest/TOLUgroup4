clc; clear; close all;

L = 30; % Lengd í mm
lambda = 0.8; % Massathettleiki g/mm^3
EI = 1.09e+10; % Pa*mm^2
% k > 0

%% Daemi 1
% plot f(x) from 0 to 10
x = linspace(0,10,100);
y = f(x);
figure
plot(x,y,'b',LineWidth=1.5)
grid on
xlabel('x')
ylabel('f(x)')
% title('f(x) = cos(x)*cosh(x) + 1')
% yticks([0 50 100]);
ylim([-10 10])
xlim([0 10])
ax = gca;
ax.XAxisLocation = 'origin';

saveas(gcf,"skil1daemi1plot.png")

%% Daemi 2
% find lowest positive root of f using bisect

% use range from 0 to 3 from plot
tol = 1e-5;

a1 = 1; b1 = 3;
k_times_L_1 = bisect(a1,b1,tol);
k1 = k_times_L_1/L;
omega_1 = omega_f(EI,lambda,k1);

a2 = 4; b2 = 6;
k_times_L_2 = bisect(a2,b2,tol);
k2 = k_times_L_2/L;
omega_2 = omega_f(EI,lambda,k2);

a3 = 7; b3 = 8;
k_times_L_3 = bisect(a3,b3,tol);
k3 = k_times_L_3/L;
omega_3 = omega_f(EI,lambda,k3);

%% Daemi 3

% Svar 4
% n: iterations/operations
% i = nr. of correct decimals + 1

% General form
% (b-a)/2^(n+1) < 10^(-i)
% 10^(-1) is chosen as a sort of degree of accuracy

% Solved for n
% n > (i + log10(b-a))/log10(2) - 1

i = 5; % We want 4 correct decimals so we chose i = nr. of correct decimals + 1
n = @(i,a,b) (i + log10(b-a))/log10(2) - 1;

% Svar 4
% n: iterations/operations
% i = nr. of correct decimals + 1

% General form
% (b-a)/2^(n+1) < 0.5*10^(-i)
% 0.5 is chosen as a sort of degree of accuracy

% Solved for n
% n > (i + log10(b-a) + log10(2))/log10(2) - 1

% i = 4; % We want 4 correct decimals so we chose i = nr. of correct decimals + 1
% n = @(i,a,b) (i + log10(b-a) + log10(2))/log10(2) - 1;

nr_of_op1 = n(i,a1,b1)
nr_of_op1 = ceil(nr_of_op1); % to round up
disp("nr. of operations for root 1 (math) = " + nr_of_op1)

nr_of_op2 = n(i,a2,b2);
nr_of_op2 = ceil(nr_of_op2); % to round up
disp("nr. of operations for root 2 (math) = " + nr_of_op2)

nr_of_op3 = n(i,a3,b3);
nr_of_op3 = ceil(nr_of_op3); % to round up
disp("nr. of operations for root 3 (math) = " + nr_of_op3)

%% Daemi 4

x0 = 4; tol = 1e-4;
k_times_L_2 = bisect(a2,b2,tol);
k2 = k_times_L_2/L;
% omega_2 = omega_f(EI,lambda,k2);

%% Daemi 5


