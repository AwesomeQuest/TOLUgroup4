function [tvec, yvec] = RKsolver(y0,n,T)
    y = zeros(length(y0),n);
    y(:,1)=y0;
    h=T/n;
    t=0;
    yvec=[y0];tvec=[0];
    for i=1:n
        t=t+h;
        k1=f(t,y(:,i));
        k2=f(t+h/2,y(:,i)+h/2.*k1);
        k3=f(t+h/2,y(:,i)+h/2.*k2);
        k4=f(t+h,y(:,i)+h.*k3);
        y(:,i+1)=y(:,i)+h/6*(k1+2.*k2+2.*k3+k4);
        tvec(i+1)=t;
        yvec(:,i+1)=y(:,i+1);
    end
end

function z=f(t,y)
    % fallið í diffurjöfnu
    alpha = 0.5; beta = 0.01;
    gamma = 0.005; delta = 0.2;

    % Differential equations
    z = [
        alpha * y(1) - beta * y(1) * y(2); % dy1/dt
        gamma * y(1) * y(2) - delta * y(2) % dy2/dt
    ];
end
