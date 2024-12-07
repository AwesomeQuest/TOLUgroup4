function [tvec, yvec] = RKsolver3(y0,n,T,b_1)
    arguments
        y0; n; T; b_1 = 3;
    end
    y = zeros(length(y0),n);
    y(:,1)=y0;
    h=T/n;
    t=0;
    yvec=[y0];tvec=[0];
    for i=1:n
        t=t+h;
        k1=ht2(t,y(:,i),b_1);
        k2=ht2(t+h/2,y(:,i)+h/2.*k1,b_1);
        k3=ht2(t+h/2,y(:,i)+h/2.*k2,b_1);
        k4=ht2(t+h,y(:,i)+h.*k3,b_1);
        y(:,i+1)=y(:,i)+h/6*(k1+2.*k2+2.*k3+k4);
        tvec(i+1)=t;
        yvec(:,i+1)=y(:,i+1);
    end
end