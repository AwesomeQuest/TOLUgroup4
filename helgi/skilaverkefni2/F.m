% function y = F(x)
%     mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;
%     T = (2:10) .*10 + 273.15;
%     T = T';
%     y = x(1) + x(2)./T + x(3)./(T.^2) - log(mu);
% end

function y = F(x)
    mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;
    T = (2:10) .*10 + 273.15;
    T = T';
    y = exp(x(1) + x(2)./T + x(3)./(T.^2)) - mu;
end