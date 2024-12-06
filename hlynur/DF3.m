function v=DF3(v0, b1)
    x=v0(1);y=v0(2);z=v0(3);
    %t=linspace(0,T,n);
    a1 = 5; a2 = 0.1; b2 = 2; d1 = 0.4; d2 = 0.01;
    dx= x*(1-x)- (a1*x*y)/(1+b1*x);
    dy= (a1*x*y)/(1+b1*x) -d1*y-(a2*y*z)/(1+b2*y);
    dz= (a2*y*z)/(1+b2*y)-d2*z;
    v= [dx; dy ;dz];
end



    