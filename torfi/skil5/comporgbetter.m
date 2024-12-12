close all; clear all;

N = 9;
T = 30;
M = 100;

btrvals = zeros(N,1);
orgvals = zeros(N,1);

[ogxs,~,~] = iterdiffv1(T,20000,M);
trorg = max(ogxs(:,end));
[trxs,~,~] = iterdiffv2(T,20000,M);
trtru = max(trxs(:,end));


for i = 1:N
    [ogxs,~,~] = iterdiffv1(T,100*2^i,M);
    orgvals(i) = max(ogxs(:,end));
    [trxs,~,~] = iterdiffv2(T,100*2^i,M);
    btrvals(i) = max(trxs(:,end));
end

xs = log2(100*2 .^(2:N));
ys1 = log2(abs(orgvals(1:end-1) - trtru) ./ (orgvals(2:end)));
ys2 = log2(abs(btrvals(1:end-1) - trtru) ./ (btrvals(2:end)));

plot(xs,ys1)
hold on
plot(xs,ys2)


fit1 = polyfit(log2(100*2 .^(2:N)), log2(abs(orgvals(1:end-1) - trtru) ./ (orgvals(2:end))), 1)
fit2 = polyfit(log2(100*2 .^(2:N)), log2(abs(btrvals(1:end-1) - trtru)./ (btrvals(2:end))), 1)

legend("Orginal: slope " + num2str(fit1(1),3) , "Better: slope " + num2str(fit2(1),4) )
xlabel("log2(N)")
ylabel("log2(\epsilon)")