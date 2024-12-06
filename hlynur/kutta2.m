function y=kutta2(y0,T,n,w)
y(:,1)=y0;h=T/n;t=0;
for i=1:n
    k1=DF2(y(:,i),w);
    k2=DF2(y(:,i)+h/2*k1,w);
    k3=DF2(y(:,i)+h/2*k2,w);
    k4=DF2(y(:,i)+h*k3,w);
    y(:,i+1)=y(:,i)+h/6*(k1+2*k2+2*k3+k4);
end
end
