clc; clear;
T = 40;

N = 200;
M = N;
k = T/N;

ts = 0:k:T;
bound = [
    ((0:(20/T*N)))*k/20, ...
    ones(1,floor(30/T*N- 20/T*N)), ...
    % 2.5 - ((30/T*N:50/T*N))*k/20, ...
    % zeros(1,floor((T-50)/T*N))
];

[ts(1:length(bound))' bound']
