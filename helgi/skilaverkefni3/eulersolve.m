function [tvec, yvec] = eulersolve(y0,n,T)
    h=T/n;
    y=y0;t=0;
    yvec=[y0];tvec=[0];
%     yvec=[y0];tvec=zeros(size(y0));
    for i=1:n
        y=y+eulerstep(h,t,y);
        t=t+h;
        tvec(i+1)=t;
        yvec(:,i+1)=y;
    end
end

function w=eulerstep(h,t,y)
    % miðar við f að neðan
    w = h * f(t,y);
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
