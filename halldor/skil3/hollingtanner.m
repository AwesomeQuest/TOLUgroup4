function dUdt = hollingtanner(t, U)
    % Constants
    a1 = 5; b1 = 3; a2 = 0.1; b2 = 2; d1 = 0.4; d2 = 0.01;
    X = U(1); Y = U(2); Z = U(3);
    
    % Equations
    dXdt = X * (1 - X) - (a1 * X / (1 + b1 * X)) * Y;
    dYdt = (a1 * X / (1 + b1 * X)) * Y - d1 * Y - (a2 * Y / (1 + b2 * Y)) * Z;
    dZdt = (a2 * Y / (1 + b2 * Y)) * Z - d2 * Z;
    
    % Output derivative
    dUdt = [dXdt; dYdt; dZdt];
end
