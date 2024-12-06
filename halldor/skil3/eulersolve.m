function [tvec, yvec]=eulersolve(y0,n,T)
h=T/n;
y=y0;t=0;
yvec=[y0];tvec=[0];
for i=1:n
    y=y+eulerstep(h,t,y);
    t=t+h;
    tvec(i+1)=t;
    yvec(:,i+1)=y;
end
end
function w=eulerstep(h,t,y)
w = h * lv(t, y);
end
