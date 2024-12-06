function y = jac2(x)
    T = (2:10) .*10 + 273.15;
    T = T';
    d1 = exp(x(1)./T + x(2).*T + x(3).*(T.^2)) .* 1./T;
    d2 = exp(x(1)./T + x(2).*T + x(3).*(T.^2)) .* T;
    d3 = exp(x(1)./T + x(2).*T + x(3).*(T.^2)) .* T.^2;
    y = [d1, d2, d3];
end