function y=kutta3(y0,T,n,b1)
y(:,1)=y0;h=T/n;t=0;
for i=1:n
    k1=DF3(y(:,i),b1);
    k2=DF3(y(:,i)+h/2.*k1,b1);
    k3=DF3(y(:,i)+h/2.*k2,b1);
    k4=DF3(y(:,i)+h.*k3,b1);
    y(:,i+1)=y(:,i)+h/6.*(k1+2*k2+2*k3+k4);
end
end
