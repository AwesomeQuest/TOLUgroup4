% function y = jac(x)
%     T = (2:10) .*10 + 273.15;
%     T = T';
%     % Sama fylki og A í dæmi 1
%     % Dálkur 1: ones
%     % Dálkur 2: 1/T
%     % Dálkur 3: 1/(T^2)
%     d1 = ones(size(T));
%     d2 = T.^(-1);
%     d3 = T.^(-2);
%     y = [d1, d2, d3];
% end

function y = jac(x)
    T = (2:10) .*10 + 273.15;
    T = T';
    d1 = exp(x(1) + x(2)./T + x(3)./(T.^2));
    d2 = exp(x(1) + x(2)./T + x(3)./(T.^2)) .* 1./T;
    d3 = exp(x(1) + x(2)./T + x(3)./(T.^2)) .* 1./(T.^2);
    y = [d1, d2, d3];
end