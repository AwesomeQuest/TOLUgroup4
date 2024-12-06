function x=gaussnewton2(x0,tol)
%til að besta lausn á F(x)=0
%forritið gerir ráð fyrir að F(x) sé skilgreint
%og sömuleiðis Jacobi fylki þess jac(x)
% oldx=x0; x = oldx - jac(x0)\(F(x0));
% counter = 0;
x=x0; oldx=x+2*tol;
while norm(x-oldx)>tol % & counter < 10000
%     disp(x)
%     disp(oldx)
%     counter = counter + 1
    oldx=x;
    J=jac2(x);
    s=J\F2(x);
%     s=(J'*J) \ (J'*F2(x)); % Otharfa lina thar sem MATLAB er lowkey goated
    x=x-s;
end
% disp("counter = " + counter)
end