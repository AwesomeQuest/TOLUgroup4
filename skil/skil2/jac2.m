function y=jac2(x)
    T = [20, 30, 40, 50, 60, 70, 80, 90, 100]+273.16;
    A=[1./T' T' (T.^2)'];
    T=T'
    exp(x(1)./T + x(2)*T + x(3)*(T.^2))
    y = A.* exp(x(1)./T + x(2)*T + x(3)*(T.^2));

end