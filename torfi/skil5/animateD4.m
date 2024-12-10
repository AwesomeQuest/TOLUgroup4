clc;close all;
T = 30;

N = 100;
M = 100;


[us1,ts,xs1] = iterdiffv1(T,N,M);
[us2,ts,xs2] = iterdiffv2(T,N,M);
[us3,ts,xs3] = iterdiffv1(T,10*N,M);

plot(xs1,us1(:,round(N/6*6)))

% for i = 1:N+1
%     plot(xs1,us1(:,i))
%     hold on
%     plot(xs2,us2(:,i))
%     plot(xs3,us3(:,10*i))
%     legend("original", "modified", "many iterations original")
%     ylim([0,1])
%     pause(0.1)
%     hold off
% end
% close all;