function y=mu2(x)
beta = [-8.83589008226628, -909.852431946823, 432616.693785727];
T0_2 = 80 + 273.15;
T1_2 = 20 + 273.15;
L = 0.1;
T = T0_2 + (T1_2 - T0_2).*x./L;
y = exp(beta(1) + beta(2)./T + beta(3)./T.^2);
end