function [us,ts,xs] = iterdiffd8(T,N,M)
D = 0.01; v = 0.1; L = 5; h = L/M; k = T/N;
alpha = k*D/h^2; beta = v*k/2/h;
ts = 0:k:T; xs = 0:h:L;

A = sparse( ...
    [2:M, 2:M, 2:M],[1:M-1, 2:M, 3:M+1], ...
    [repelem(-beta-alpha,M-1), ...
     repelem(1+2*alpha,  M-1), ...
     repelem(beta-alpha, M-1), ...
    ],M+1,M+1 ...
);
A(1,1) = 1; A(end,end) = 3; A(end,end-1) = -4; A(end,end-2) = 1;

A2 = sparse( ...
    [2:M, 2:M, 2:M],[1:M-1, 2:M, 3:M+1],...
    [repelem(-2*beta-2*alpha,M-1), ...
     repelem(3+4*alpha,      M-1), ...
     repelem(2*beta-2*alpha, M-1), ...
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