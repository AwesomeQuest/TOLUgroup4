clc;close all;
N = 100;
M = N;
T = 120;

[us1,ts,xs1] = iterdiffd8(T,N,M);

[M,I] = max(us1,[], 1);

yyaxis right
plot(ts,M)
ylim([-0.02,1.01])
hold on
ylabel("C_{max}(t)")
yyaxis left
plot(ts,xs1(I))
ylim([-0.1,5.05])
legend("C(t)", "x(t)")
xlabel("t")
ylabel("x_{max}(t)")
hold off