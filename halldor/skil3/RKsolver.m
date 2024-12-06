function y=RKsolver(y0,T,n)
y(1,:)=y0';
h=T/n;
t=0;
for i=1:n
    t=t+h;
    k1=lv(t,y(i,:)');
    k2=lv(t+h/2,y(i,:)'+h/2*k1);
    k3=lv(t+h/2,y(i,:)'+h/2*k2);
    k4=lv(t+h,y(i,:)'+h*k3);
    y(i+1,:)=y(i,:)+h/6*(k1'+2*k2'+2*k3'+k4');
end
end

