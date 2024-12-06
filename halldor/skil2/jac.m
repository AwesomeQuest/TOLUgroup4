function y=jac(x)
T = (2:10) .*10 + 273.15;
T = T';
d1 = ones(size(T)).*exp(x(1) + x(2) ./T + x(3)./T.^2);
d2 = T .^-1.*exp(x(1) + x(2) ./T + x(3)./T.^2);
d3 = T.^-2.*exp(x(1) + x(2) ./T + x(3)./T.^2);
y = [d1,d2,d3];
end
 