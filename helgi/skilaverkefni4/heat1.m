function w=heat1(T,m,n)
    D=4;L=10;
    h=L/m;k=T/n;sigma=(D*k)/h^2;
    i=[2:m , 2:m , 2:m]';
    j=[1:m-1 , 2:m , 3:m+1]';
    values=[-sigma*ones(m-1,1);(1+2*sigma)*ones(1,1);-sigma*ones(m-1,1)];
    A=sparse(i,j,values);
    A(1,1)=1;A(m+1,m+1)=1;
    
    b=20*ones(m+1,1);b(1)=50;
    w(:,1)=20*ones(m+1,1); %initial value
    for j=1:n
        w(:,j+1)=A\b;
        b=[50;w(2:m,j+1);20];
    end
end