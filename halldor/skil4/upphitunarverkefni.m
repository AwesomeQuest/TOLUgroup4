% Upphitunarverkefni vika 3
% part 1
clc;clear all;
n = 10; 
h = 1 / n;
num_points = n - 1; 

A = zeros(num_points, num_points);

for i = 1:num_points
    A(i, i) = -(2 + h^2); 
    if i > 1
        A(i, i-1) = 1; 
    end
    if i < num_points
        A(i, i+1) = 1; 
    end
end

disp('Matrix A:');
disp(A);

%% part 2 

b = zeros(num_points, 1);
y0 = -1; 
yn = 2;
b(1) = -y0;
b(end) = -yn;

y_inner = A\b; % y = A^(-1)*b

y = [y0; y_inner; yn];

x = linspace(0, 1, n+1);
exact = 1.0074 * exp(x) - 2.0074 * exp(-x);

hold on
plot(x, y, 'o-', 'LineWidth', 1.5);
plot(x, exact, 'r--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('Solution with Dirichlet Boundary Conditions');

%% part 3
clc;clear all; 
n = 10; 
h = 1 / n;
A = zeros(n+1, n+1);
b = zeros(n+1, 1);

for i = 2:n
    A(i, i) = -(2 + h^2); % Main diagonal
    A(i, i-1) = 1;        % Subdiagonal
    A(i, i+1) = 1;        % Superdiagonal
end

A(1, 1:3) = [-3 / (2 * h), 4 / (2 * h), -1 / (2 * h)];
A(end, end-2:end) = [1 / (2 * h), -4 / (2 * h), 3 / (2 * h)];

b(1) = 2 * h; 
b(end) = -6 * h; 

y = A \ b;


x = linspace(0, 1, n+1);
exact = -1.43289 * exp(x) - 2.43289 * exp(-x);

figure;
plot(x, y, 'o-', 'LineWidth', 1.5); hold on;
plot(x, exact, 'r--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('Solution with Neumann Boundary Conditions');
legend('Numerical Solution', 'Exact Solution', 'Location', 'Best');

%% part 4
