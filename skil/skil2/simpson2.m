function I=simpson2(a,b,n)
h=(b-a)/n;
x=a;y=(1/mu2(a));
for i=1:n-1
    x=x+h;
    y=y+2*(1/mu2(x));
end
x=a-h/2;
for i=1:n
    x=x+h;
    y=y+4*(1/mu2(x));
end
y=y+(1/mu2(b));
I=h/6*y;
end
