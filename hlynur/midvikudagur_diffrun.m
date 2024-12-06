f=@(x) log(x); 
x0=2;
raungildi=1/2;
h=10^(-5);

for i =1:10
    f(x0+h)
    f(x0-h)
    nalgun=(f(x0+h)-f(x0-h))/(2*h);
    error(i)=abs(nalgun-raungildi);
    h=h/10;
end

