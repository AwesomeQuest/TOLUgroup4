clc; clear;
T = 100;

N = 200;
M = N;
k = T/N;

ts = 0:k:T;
bound = [
    ((0:(20/T*N)))*k/20, ...
    ones(1,30/T*N- 20/T*N), ...
    2.5 - ((30/T*N:50/T*N-1))*k/20, ...
    zeros(1,(T-50)/T*N)
];

[ts' bound']