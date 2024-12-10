function [tvec, yvec] = RKsolver2(y0,n,T,w)
    arguments
        y0; n; T; w = 0.3;
    end
    y = zeros(length(y0),n);
    y(:,1)=y0;
    h=T/n;
    t=0;
    yvec=[y0];tvec=[0];
    for i=1:n
        t=t+h;
        k1=ht(t,y(:,i),w);
        k2=ht(t+h/2,y(:,i)+h/2.*k1,w);
        k3=ht(t+h/2,y(:,i)+h/2.*k2,w);
        k4=ht(t+h,y(:,i)+h.*k3,w);
        y(:,i+1)=y(:,i)+h/6*(k1+2.*k2+2.*k3+k4);
        tvec(i+1)=t;
        yvec(:,i+1)=y(:,i+1);
    end
end