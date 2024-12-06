function dydt = lv(t, y)
alpha=0.5;beta=0.01;gamma=0.005;delta=0.2;
dydt = [alpha.*y(1) - beta.*y(1).*y(2);
        gamma.*y(1).*y(2) - delta.*y(2)];
end