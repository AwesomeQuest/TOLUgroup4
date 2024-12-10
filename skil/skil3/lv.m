function z = lv(t,y)
    % fallið í diffurjöfnu
    alpha = 0.5; beta = 0.01;
    gamma = 0.005; delta = 0.2;

    z = [
        alpha * y(1) - beta * y(1) * y(2); % dy1/dt
        gamma * y(1) * y(2) - delta * y(2) % dy2/dt
    ];
end