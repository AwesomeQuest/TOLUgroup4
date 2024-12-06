function [tvec, yvec] = RK4solveht(y0, T, n)
    % RK4solve með holling tanner
    h = T / n;
    t = 0;

    yvec = zeros(n+1, length(y0)); 
    tvec = linspace(t, T, n+1); 

    yvec(1, :) = y0;

    for i = 1:n
        y = yvec(i, :)';

        k1 = h * ht(tvec(i), y);
        k2 = h * ht(tvec(i) + h / 2, y + k1 / 2);
        k3 = h * ht(tvec(i) + h / 2, y + k2 / 2);
        k4 = h * ht(tvec(i) + h, y + k3);

        y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

        yvec(i+1, :) = y_next';
    end
end