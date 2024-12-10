function [tvec, yvec] = RK4solveCLV(y0,n,T)
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

function z=f(t,x)
    % fallið í diffurjöfnu
    alpha = [
        0.238238,  0.195195,  0.219219;
        0.348348,  0.244244,  0.145145;
        0.232232,  0.252252,  0.200200
    ];
    r = ones(3,1);

    % Differential equations
    z = r .* x .* ( 1 - alpha * x);
end
