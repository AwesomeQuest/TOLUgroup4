function y=mu_func2(T,beta)
    y= exp(beta(1)./T' + beta(2)*T' + beta(3)*(T.^2)');
end
