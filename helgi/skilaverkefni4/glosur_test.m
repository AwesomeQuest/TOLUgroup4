% L=10; 
clc; clear; close all;
T=10; m=100; n=100;
% h=L/m; k=T/n;
w=heat1(T,m,n);
%plot
x=(0:m)*h; t=(0:n)*k;
mesh(x,t,w');

tracerod = animatedline('Color','b','LineWidth',1);
for j = 1:n+1
    clearpoints(tracerod)
    for i=1:m+1
        addpoints(tracerod,x(i),w(i,j));
    end
    %pause(0.1) use to slow down animation
    drawnow
end

% D=4;L=10;
% h=L/m;k=T/n;sigma=(D*k)/h^2;
% i=[2:m , 2:m , 2:m]';
% j=[1:m-1 , 2:m , 3:m+1]';
% values=[-sigma*ones(m-1,1);(1+2*sigma)*ones(1,1);-sigma*ones(m-1,1)];
% A=sparse(i,j,values);
% A(1,1)=1;A(m+1,m+1)=1;
