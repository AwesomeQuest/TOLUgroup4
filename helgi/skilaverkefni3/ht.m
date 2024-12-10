function z = ht(t,y,w)
    arguments
        t; y; w = 0.3;
    end
    r_1 = 1; r_2 = 0.1; k = 7;
    d = 1; j = 1;
    z = [
        r_1*y(1)*(1 - (y(1)/k)) - (w*y(1))/(d + y(1))*y(2);
        r_2*y(2)*(1 - (j*y(2))/y(1))
    ];
end