clc; clear;
T = 50;

N = 220;
M = N;
k = T/N;

ts = 0:k:T;
bound = arrayfun(@(j) f(j,T,N), 1:N+1);


[ts(1:length(bound))' bound']


function C = f(j,T,N)
    k = T/N;
    if j <= 20/T*N
        C = (j-1)*k/20;
    elseif j >= 20/T*N && j <= 30/T*N
        C = 1;
    elseif j >= 30/T*N && j <= 50/T*N
        C = 2.5 - (j-1)*k/20;
    else
        C = 0;
    end
end