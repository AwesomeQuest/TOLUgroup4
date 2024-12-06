function x=newton(x0,tol)
%skilgreina hér fallið f(x)
%f=@(x) ...
%og f'(x)
%Df=@(x) ...
x=x0;oldx=x-100; % oldx = x0 + 2*tol
counter=0; maxiterations=30;
while abs(x-oldx)>tol
    oldx=x;
    x=x-f(x)/Df(x);
    counter = counter + 1;
    if counter > maxiterations
        disp("Newton method unsuccessful")
        break
    end
end
end
