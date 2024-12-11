clc;close all;
T = randi([1,500]);

N = 30;
M = N;


[us1,ts1,xs1] = iterdiffd8(T,N,M);
% [us1,ts,xs1] = iterdiffd8org(T,N,M);

% plot(xs1,us1(:,round(N/6*6)))

for i = 1:N+1
    plot(xs1,us1(:,i))
    text(2.5,2.5,num2str(ts(i)))
    ylim([0,1])
    pause(0.1)
end