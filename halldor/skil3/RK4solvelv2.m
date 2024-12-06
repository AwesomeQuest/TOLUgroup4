function [tvec, yvec] = RK4solvelv2(y0, n, T)
    % RK4solve með RK4solve
    h = T / n;
    t = 0;
    y = y0;
    
    tvec = zeros(1, n+1);
    yvec = zeros(length(y0), n+1);

    tvec(1) = t;
    yvec(:, 1) = y0;

    for i = 1:n
        k1 = h * lv(t, y);
        k2 = h * lv(t + h / 2, y + k1 / 2);
        k3 = h * lv(t + h / 2, y + k2 / 2);
        k4 = h * lv(t + h, y + k3);

        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

        t = t + h;

        tvec(i+1) = t;
        yvec(:, i+1) = y;
    end
end
