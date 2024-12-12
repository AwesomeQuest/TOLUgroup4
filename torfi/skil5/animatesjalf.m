% clc;close all;
N = 100;
M = 400;
T = 120;

[us1,ts1,xs1] = iterdiffsjalf(T,N,M);

% plot(xs1,us1(:,round(N/6*6)))

for i = 1:N+1
    plot(xs1,us1(:,i))
    text(2.5,0.5,num2str(ts1(i)))
    ylim([0,1])
    xlim([0,5])
    pause(0.1)
end