function z = ht2(t,y,b_1)
    arguments
        t; y; b_1 = 3;
    end
%     if nargin < 4 % Ef nr of args er minna en 4 er b_1 fastinn settur sem 3
%         b_1 = 3; % gerir mer kleift ad profa f mismunandi b_1
%     end % akvad ad nota frekar fetusinn ad ofan
    a_1 = 5; a_2 = 0.1; % b_1 = 3;
    b_2 = 2; d_1 = 0.4; d_2 = 0.01;
    z = [
        y(1)*(1 - y(1)) - (a_1*y(1)*y(2))/(1 + b_1*y(1));
        (a_1*y(1)*y(2))/(1 + b_1*y(1)) - d_1*y(2) - (a_2*y(2)*y(3))/(1 + b_2*y(2));
        (a_2*y(2)*y(3))/(1 + b_2*y(2)) - d_2*y(3)
    ];
end