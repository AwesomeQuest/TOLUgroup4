T = 30;

M = 10
N = 3

D = 0.01;
v = 0.1;
L = 5;
h = L/M;
k = T/N;

ts = 0:k:T;

alpha = k*D/h^2;
beta  = v*k/2/h;


us = zeros(M+1,N+1);
for i = 2:M
    us(i,1) = exp(-(h*(i-1) -1)^2/D);
end
us(:,1)