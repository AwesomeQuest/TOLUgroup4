%% Daemi 2
clc; clear; close all;

% Jadargildi: y(0)=-1, y(1)=2

y_nakvaem1 = @(x) 1.0074*exp(x) - 2.0074*exp(-x);

n = 10; h = 1/(n); 
alpha = 1; beta = -(2 + h^2); gamma = 1;

A1 = zeros(n+1);

A1(1,1) = 1;
A1(n+1,n+1) = 1;

for i = 1:n-1
    A1(i+1,i) = alpha;
    A1(i+1,i+1) = beta;
    A1(i+1,i+2) = gamma;
end

b1 = zeros(n+1,1);
b1(1) = -1; b1(end) = 2;

x1 = A1\b1;

t1 = linspace(0,1,n+1);
y_real1 = y_nakvaem1(t1);

hold on
plot(t1,y_real1,"b",LineWidth=1.5);
plot(t1,x1,"--r",LineWidth=1.5);

%% Daemi 3
clc; clear; close all;

% Jadargildi: y'(0)=1, y'(1)=-3

y_nakvaem2 = @(x) -1.43289*exp(x) - 2.43289*exp(-x);

n = 10; h = 1/(n);
alpha = 1; beta = -(2 + h^2); gamma = 1;
A2 = zeros(n+1);

A2(1,1) = -3; A2(1,2) = 4; A2(1,3) = -1;
A2(end,end-2) = 1; A2(end,end-1) = -4; A2(end,end) = 3;

for i = 1:n-1
    A2(i+1,i) = alpha;
    A2(i+1,i+1) = beta;
    A2(i+1,i+2) = gamma;
end

b2 = zeros(n+1,1);
b2(1) = 2*h; b2(end) = -6*h

x2 = A2\b2;

t2 = linspace(0,1,n+1);
y_real2 = y_nakvaem2(t2);

hold on
plot(t2,y_real2,"b",LineWidth=1.5);
plot(t2,x2,"--r",LineWidth=1.5);

%% Daemi 4
clc; clear; close all;

% Jadargildi: y'(0)=2*y(0), y'(1)=1-y(1)

y_nakvaem3 = @(x) 0.18394*exp(x) - 0.06131*exp(-x);

n = 100; h = 1/(n-1);
alpha = 1; beta = -(2 + h^2); gamma = 1;
A3 = zeros(n+1);

A3(1,1) = -(3+4*h); A3(1,2) = 4; A3(1,3) = -1;
A3(end,end-2) = 1; A3(end,end-1) = -4; A3(end,end) = 3+2*h;

for i = 1:n-1
    A3(i+1,i) = alpha;
    A3(i+1,i+1) = beta;
    A3(i+1,i+2) = gamma;
end

b3 = zeros(n+1,1); 
b3(end) = 2*h;

x3 = A3\b3;

t3 = linspace(0,1,n+1);
y_real3 = y_nakvaem3(t3);

hold on
plot(t3,y_real3,"b",LineWidth=1.5);
plot(t3,x3,"--r",LineWidth=1.5);

