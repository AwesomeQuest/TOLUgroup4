function y=F(u)
    mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]*1e-3;
    T = [20, 30, 40, 50, 60, 70, 80, 90, 100]+273.16;
    y= exp(u(1) +u(2)./T' + u(3)./(T.^2)')-mu';
end
