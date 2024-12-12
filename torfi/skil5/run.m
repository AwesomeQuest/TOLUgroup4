
T = 30;

N = 100;
M = 100;


[us,ts] = iterdiffv1(T,3,10)

plot(us(:,1))

function [us,ts] = iterdiffv1(T,N,M)
    D = 0.01;
    v = 0.1;
    L = 5;
    h = L/M;
    k = T/N;

    ts = 0:k:T;

    alpha = k*D/h^2;
    beta  = v*k/2/h;

    A = sparse(M+1,M+1);
    A(1,1) = 1;
    A(end,end) = 1;
    for i = 2:M
        A(i,i-1) = -beta-alpha;
        A(i,i) = 1+2*alpha;
        A(i,i+1) = beta-alpha;
    end

    us = zeros(M+1,N+1);
    for i = 2:M
        us(i,1) = exp(-(h*(i-1) -1)^2/D);
    end
    
    for i = 2:N+1
        us(:,i) = A\us(:,i-1);
    end
    return
end