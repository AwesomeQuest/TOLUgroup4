function y = mu_func1(T,beta)
    y = exp(beta(1) + beta(2)./T + beta(3)./(T.^2));
end