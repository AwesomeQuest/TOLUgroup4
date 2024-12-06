function x=gaussnewton(x0,tol)
%til að besta lausn á F(x)=0
%forritið gerir ráð fyrir að F(x) sé skilgreint
%og sömuleiðis Jacobi fylki þess jac(x)
x=x0;oldx=x+2*tol;
while norm(x-oldx)>tol
    oldx=x;
    J=jac(x);
    s=(J'*J) \ (J'*F(x));
    x=x-s;
end
end
