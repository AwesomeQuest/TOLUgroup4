function [us,ts,xs] = iterdiffsjalf2(T,N,M)
D = 0.01; L = 5; h = L/M; k = T/N;
ts = 0:k:T; xs = 0:h:L;

alpha1 = @(j) -v(j*h)*k/2/h - k*D/h^2;
beta1 = @(j) 1+2*k*D/h^2 - k*vp(j*h);
gamma1 = @(j) v(j*h)*k/2/h - k*D/h^2;

alpha2 = @(j) -2*k*(D/h^2 + v(j*h)/h/2);
beta2 = @(j) -(-3-2*k*(vp(j*h) + 2*D/h^2));
gamma2 = @(j) -2*k*(D/h^2 - v(j*h)/2/h);

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

us = zeros(M+1,N+1);
bound = arrayfun(@(j) f(j,T,N), 1:N+1);
us(1,:) = bound;
us(:,2) = A\us(:,1);
us(1,2) = bound(2);

bv = us(:,2);
bv(end) = 0;
bv(2:M)  = 4*us(2:M,2) - us(2:M,1);
for j = 3:N+1
    us(:,j) = A2\bv;
    us(1,j) = bound(j);
    bv(1) = bound(j); bv(end) = 0;
    bv(2:M) = 4*us(2:M,j) - us(2:M,j-1);
end
end

function C = f(j,T,N)
    k = T/N;
    if j <= 20/T*N
        C = (j-1)*k/20;
    elseif j >= 20/T*N && j <= 30/T*N
        C = 1;
    elseif j >= 30/T*N && j <= 50/T*N
        C = 2.5 - (j-1)*k/20;
    else
        C = 0;
    end
end

function y = v(x)
    y = (x^2+6*x+3)/600;
end
function y = vp(x)
    y = (1/600)*(6 + 2*x);
end