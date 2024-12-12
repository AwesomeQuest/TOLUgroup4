clc; clear;
T = 50;

N = 220;
M = N;
k = T/N;

ts = 0:k:T;


alpha1 = @(j) -v(j*h)*k/2/h - k*D/h^2;
beta1 = @(j) 1+2*k*D/h^2 + k*vp(j*h);
gamma1 = @(j) v(j*h)*k/2/h - k*D/h^2;

alpha2 = @(j) 2*k*(D/h^2 + vp(j*h)/h);
beta2 = @(j) (-3-2*k*(vp(j*h) + 2*D/h^2));
gamma2 = @(j) 2*k*(D/h^2 - vp(j*h)/2/h);

A = sparse( ...
    [2:M, 2:M, 2:M],[1:M-1, 2:M, 3:M+1], ...
    [arrayfun(alpha1, 1:M-1), ...
     arrayfun(beta1, 1:M-1), ...
     arrayfun(gamma1, 1:M-1), ...
    ],M+1,M+1 ...
);
A(1,1) = 1; A(end,end) = 3; A(end,end-1) = -4; A(end,end-2) = 1;

A2 = sparse( ...
    [2:M, 2:M, 2:M],[1:M-1, 2:M, 3:M+1],...
    [arrayfun(alpha2, 1:M-1), ...
     arrayfun(beta2, 1:M-1), ...
     arrayfun(gamma2, 1:M-1), ...
    ], M+1,M+1 ...
);
A2(1,1) = 1; A2(end,end) = 3; A2(end,end-1) = -4; A2(end,end-2) = 1;


function y = v(x)
    y = (x^2+6*x+3)*sin(x^2)^2/60 + 0.00001;
end

function y = vp(x)
    y = (1/60)*(6 + 2*x)*(sin(x^2)^2) + (1/30)*x*(3 + 6*x + x^2)*sin(2*(x^2));
end
