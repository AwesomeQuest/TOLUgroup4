function y = T(x)
    T0 = 20 + 273.15;
    T1 = 80 + 273.15;
    L = 0.1;
    y = T0 + x*(T1-T0)/L;
end