function y=jac2(x)
T = (2:10) .*10 + 273.15;
T = T';
d1 = 1./T.*exp(x(1)./T + x(2) .*T + x(3).*T.^2);
d2 = T .*exp(x(1)./T + x(2) .*T + x(3).*T.^2);
d3 = T.^2.*exp(x(1)./T + x(2) .*T + x(3).*T.^2);
y = [d1,d2,d3];
end