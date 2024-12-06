clc
T = [20, 30, 40, 50, 60, 70, 80, 90, 100]+273.16;
mu = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]*1e-3;

A = [ones(size(T')) 1./T' 1./(T.^2)']

x = A \ log(mu')

% Verkefni 2
tol=1e-6;
s=gaussnewton([0,0,0]',tol)
first_error=norm(A*x-log(mu'))
second_error=norm(F(s))

% Verkefni 3
B=[ 1./T' T' (T.^2)'];
x_3_1= B \log(mu')
norm(B*x_3_1-log(mu'))
x_3_2= gaussnewton2([0,0,0]',tol)
norm(mu_func2(T,x_3_2)-mu')
norm(mu_func2(T,x_3_1)-mu')
%%
% Verkefni 4

v0=0.01;
L=0.1;
n=2;
y=linspace(0,L,100);
vel=zeros(size(y));
for i =1:length(y)
    current_velocity=velocity(y(i),L,v0,n);
    vel(i)= current_velocity;
end
plot(y,vel)



%% Dæmi 5
v0=0.01;
L=0.1;
n = [2, 4, 8, 16, 32, 64, 128, 256];
y=L/2;
velocity_list=zeros(size(n));
for i=1:length(n)
    velocity_list(i)=velocity(y,L,v0,n(i));
end
real_velocity=velocity(y,L,v0,2^14)
%reiknum error á milli falla (mun á hraða)
error=real_velocity-velocity_list
factor_munur= error(1:end-1)./error(2:end)
log2_munur=log2(factor_munur)


%% Dæmi 6

V0 = 0.01; tol = 0.5e-10; a = 0; L = 0.1; y2 = L/2;

nefnari = adapquad(a,y2,tol, "numerator");
teljari = adapquad(a,L,tol, "denominator");

V = V0 * (nefnari/teljari)
adapquad();


%% Dæmi 7


ys = 0:0.005:L;
vel=zeros(size(ys));
n = 256;
for i =1:length(ys)
    current_velocity=velocity(ys(i),L,v0,n);
    vel(i)= current_velocity;
end

plot(ys,10*vel)
hold on;
scatter_data = zeros(size(ys));
scatter_plot = scatter(zeros(size(ys)), ys, 'r');

scatter_plot.XDataSource = 'scatter_data';
fps = 30;
tend = 10;
dt = tend / (fps * tend); 
for t = 0:dt:tend
    scatter_data = vel * t;
    refreshdata
    drawnow
    pause(1 / fps); 
end



%% Dæmi 7


ys = 0:0.005:L;
vel=zeros(size(ys));
n = 256;
for i =1:length(ys)
    current_velocity=velocity(ys(i),L,v0,n);
    vel(i)= current_velocity;
end

plot(ys,10*vel)
hold on;
scatter_data = zeros(size(ys));
scatter_plot = scatter(zeros(size(ys)), ys, 'r');

scatter_plot.XDataSource = 'scatter_data';
fps = 30;
tend = 10;
dt = tend / (fps * tend); 
for t = 0:dt:tend
    scatter_data = vel * t;
    refreshdata
    drawnow
    pause(1 / fps); 
end


%% Dæmi 8


ys = 0:0.005:L;
vel=zeros(size(ys));
n = 256;
for i =1:length(ys)
    current_velocity=velocity2(ys(i),L,v0,n);
    vel(i)= current_velocity;
end

plot(ys,10*vel)
hold on;
scatter_data = zeros(size(ys));
scatter_plot = scatter(zeros(size(ys)), ys, 'r');

scatter_plot.XDataSource = 'scatter_data';
fps = 30;
tend = 10;
dt = tend / (fps * tend); 
for t = 0:dt:tend
    scatter_data = vel * t;
    refreshdata
    drawnow
    pause(1 / fps); 
end