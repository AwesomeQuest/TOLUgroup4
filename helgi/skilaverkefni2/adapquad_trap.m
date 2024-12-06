function x=adapquad_trap(a,b,tol)
%skilgreina hér fallið g(t) sem er verið að heilda
c=(a+b)/2;
sab=trap(a,b);sac=trap(a,c);scb=trap(c,b);
if abs(sab-sac-scb)<3*tol %ath. 3 passar við trapísuaðferð
    x=sac+scb;
else
    x=adapquad_trap(a,c,tol/2)+adapquad_trap(c,b,tol/2);
end
end

function s=trap(a,b)
s=(g(a)+g(b))*(b-a)/2;
end
