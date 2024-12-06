function x=adapquad(a,b,tol)
    c=(a+b)/2;
    sab=simpson(a,b,1);sac=simpson(a,c,1);scb=simpson(c,b,1);
    if abs(sab-sac-scb)<10*tol 
        x=sac+scb;
    else
        x=adapquad(a,c,tol/2)+adapquad(c,b,tol/2);
    end
end



