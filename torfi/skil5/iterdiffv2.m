function [us,ts,xs] = iterdiffv2(T,N,M)
    D = 0.01;
    v = 0.1;
    L = 5;
    h = L/M;
    k = T/N;

    ts = 0:k:T;
    xs = 0:h:L;

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

    A2 = sparse(M+1,M+1);
    A2(1,1) = 1;
    A2(end,end) = 1;
    for i = 2:M
        A2(i,i-1) = -2*beta-2*alpha;
        A2(i,i) = 3+4*alpha;
        A2(i,i+1) = 2*beta-2*alpha;
    end

    us = zeros(M+1,N+1);
    for i = 2:M
        us(i,1) = exp(-(h*(i-1) -1)^2/D);
    end
    us(:,2) = A\us(:,1);
    us(1,2) = 0;
    us(end,2) = 0;

    bv = us(:,2);
    bv(1) = 0;
    bv(end) = 0;
    for i = 2:M
        bv(i)  = 4*us(i,2) - us(i,1);
    end
    for j = 3:N+1
        us(:,j) = A2\bv;
        us(1,j) = 0;
        us(end,j) = 0;
        for i = 1:M+1
            bv(i)  = 4*us(i,j) - us(i,j-1);
        end
    end
    return
end