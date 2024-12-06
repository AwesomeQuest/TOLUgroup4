%% Skilaverkefni 3

%upphafsfastar
alpha=0.5; beta=0.01; gamma =0.005; delta=0.2;
axis_max=100;
[y1,y2]=meshgrid(0:5:axis_max,0:5:axis_max);
dy1=Dy1(y1,y2);
dy2=Dy2(y1,y2);
normalized=sqrt(y1.^2+y2.^2);
diff_norm=sqrt(dy1.^2+dy2.^2);
y1_normalized= y1 ./normalized;
y2_normalized= y2 ./ normalized;
dy1_norm= dy1 ./diff_norm;
dy2_norm= dy2 ./diff_norm;

quiver(y1,y2,dy1_norm,dy2_norm);
axis([0 axis_max 0 axis_max]);
grid on


%% Dæmi 2
y1_0=30;y2_0=3;
y0=[y1_0;y2_0];
T=50;n=100;emat=eulersolve(y0,T,n);
plot(emat(1,:),emat(2,:),"g"),
hold on 
plot(emat(1,:),emat(3,:),"r")
hold off
%% Dæmi 3
quiver(y1,y2,dy1,dy2)
hold on
plot(emat(1,:),emat(2,:))
hold off
%% Dæmi 4
clc
close all
y1_0=40;y2_0=9;
y0=[y1_0;y2_0];
T=50;
tspan=[0 T];
%n=[100 200 400 800];
n=[100 200 400 800 1600,3200,6400,12800,25600,51200,102400];
kutta_val=kutta(y0,T,n(end)*2);
error=zeros(2,length(n));
for i=1:length(n)
    rec_emat=eulersolve(y0,T,n(i));
    %norm_vectors(:,i)=[interpolated_norm(rec_emat(2,:), kutta_val(1,:)); interpolated_norm(rec_emat(3,:),kutta_val(2,:))];
    error(:,i)= [abs(rec_emat(2,end)-kutta_val(1,end)) ;abs(rec_emat(3,end)-kutta_val(2,end))];
end

fractional_error_euler=error(:,1:end-1)./error(:,2:end);


%test_norm=norm_vectors(:,1:end-1)./norm_vectors(:,2:end);

%% Dæmi 5
T=50;n=800;
y1_0=30;y2_0=5;
y0=[y1_0;y2_0];
rk=kutta(y0,T,n);

real_val=rk
n=[100 200 400 800];
error_euler=zeros(2,length(n));
for i=1:length(n)
    rec_emat=eulersolve(y0,T,n(i));
    %error(:,i)= [abs(rec_emat(2,end)-y(1,end)) ;abs(rec_emat(3,end)-y(2,end))];
    error_euler(:,i)=[abs(rec_emat(2,end)-real_val(1,end)) ; abs(rec_emat(3,end)-real_val(2,end))];
end
disp(error_euler)
plot(n,error_euler)

% plot(rec_emat(1,:), rec_emat(2:end,:), "r", LineWidth=1)
% hold on
% plot(rec_emat(1,:),rk, "g", LineWidth=1 )
% legend("Euler y1", "Euler y2", "RK y1", "RK y2")

%% Dæmi 5 aftur
y0_1 = [10; 10]; % 1
y0_2 = [4; 4]; % 2
y0_3 = [2; 2]; % 3
T = 20;
n = 500;
t=linspace(0,T,n+1);
rk_1=kutta(y0_1,T,n);
rk_2=kutta(y0_2,T,n);
rk_3=kutta(y0_3,T,n);
hold on 
plot(t, rk_1(1, :), 'b', t, rk_1(2, :), '--b', 'LineWidth', 1.5);
plot(t, rk_2(1, :), 'g', t, rk_2(2, :), '--g', 'LineWidth', 1.5); % dæmi 2
plot(t, rk_3(1, :), 'r', t, rk_3(2, :), '--r', 'LineWidth', 1.5); % dæmi 3
legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2','Location', 'northwest');
xlabel('Tími');
ylabel('Stofn');

title("Lotka-Volterra líkan leyst með aðferð Runge-Kutta (RK4)")
hold off

% plot(rk_1(1,:),rk_1(2,:),'b',LineWidth=1.1)
% hold on
% plot(rk_2(1,:),rk_2(2,:),'g',LineWidth=1.1)
% plot(rk_3(1,:),rk_3(2,:),'r',LineWidth=1.1)
% legend("y_0=[10 10]", "y_0=[4 4]", "y_0=[2 2]")
% xlabel('Bráð');
% ylabel('Rándýr');
% title("Phase plot Lotka-Volterra líkan leyst með aðferð Runge-Kutta (RK4)");
%% Dæmi 6
%n=[3200,6400,12800,256000,512000];
n=[100,200,400,800,1600,3200,6400,12800];
%n=[1600,3200,6400,12800];
T=50 ;
y1_0=40;y2_0=9;
%y1_0=50;y2_0=50;
y0=[y1_0;y2_0];
real_n=1e6; %1.000.000
real_value=kutta(y0,T,real_n);
error_rk=zeros(2,length(n));
for i=1:length(n)
    rk=kutta(y0,T,n(i));
    error_rk(:,i)=abs(rk(:,end)-real_value(:,end));
end

change_in_error=error_rk(:,1:end-1)./ error_rk(:,2:end)
loglog(n,error_rk)
%plot(n,error_rk)
%fprintf("%1.0f  ",change_in_error)

%% Hluti 2 
% dæmi 1
y0=[10;2];
n=3200;T=50; w=0.3;
rk=kutta2(y0,T,n,w);
t=linspace(0,T,n+1);

plot(t, rk(1, :), 'b', t, rk(2, :), '--b', 'LineWidth', 1.5);
xlabel('Tími');
ylabel('Stofn');
legend('Bráð y_0= 10', 'Rándýr y_0= 2');

%% 
% dæmi 2
y=[6;4]
y1=y(1);
y2=y(2);
r1 = 1; r2 = 0.1; k = 7; d = 1; j = 1; w = 0.3;
dy1=r1*y1*(1-y1/k)-w*y1*y2/(d+y1)
dy2=r2*y2*(1-((j*y2)/y1))

%% dæmi 2 aftur
%% dæmi 3
y0_1=[10;10];
y0_2=[4;4];
y0_3=[2;2];
n=600;T=42; w=0.3;
rk_1=kutta2(y0_1,T,n,w);
rk_2=kutta2(y0_2,T,n,w);
rk_3=kutta2(y0_3,T,n,w);
t=linspace(0,T,n+1);
hold on

plot(t, rk_1(1, :), 'b', t, rk_1(2, :), '--b', 'LineWidth', 1.5);
plot(t, rk_2(1, :), 'g', t, rk_2(2, :), '--g', 'LineWidth', 1.5); % dæmi 2
plot(t, rk_3(1, :), 'r', t, rk_3(2, :), '--r', 'LineWidth', 1.5); % dæmi 3

legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2');
xlabel('Tími');
ylabel('Stofn');
title("Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")
hold off
%% dæmi 3
y0_1=[10;10];
y0_2=[4;4];
y0_3=[2;2];
n=600;T=120; w=1;
rk_1=kutta2(y0_1,T,n,w);
rk_2=kutta2(y0_2,T,n,w);
rk_3=kutta2(y0_3,T,n,w);
t=linspace(0,T,n+1);
hold on

plot(t, rk_1(1, :), 'b', t, rk_1(2, :), '--b', 'LineWidth', 1.5);
plot(t, rk_2(1, :), 'g', t, rk_2(2, :), '--g', 'LineWidth', 1.5); % dæmi 2
plot(t, rk_3(1, :), 'r', t, rk_3(2, :), '--r', 'LineWidth', 1.5); % dæmi 3

legend('Bráð (tilfelli 1) y_0= 10', 'Rándýr (tilfelli 1) y_0= 10', 'Bráð (tilfelli 2) y_0= 4', 'Rándýr (tilfelli 2) y_0= 4', 'Bráð (tilfelli 3) y_0= 2', 'Rándýr (tilfelli 3) y_0= 2');
xlabel('Tími');
ylabel('Stofn');
title("Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")
hold off
% 
% plot(rk_1(1,:),rk_1(2,:),'b',LineWidth=1.1)
% hold on
% plot(rk_2(1,:),rk_2(2,:),'g',LineWidth=1.1)
% plot(rk_3(1,:),rk_3(2,:),'r',LineWidth=1.1)
% legend("y_0=[10 10]", "y_0=[4 4]", "y_0=[2 2]")
% xlabel("Bráð");
% ylabel('Rándýr');
% title("Phase plot Holling-Tanner líkan leyst með aðferð Runge-Kutta (RK4)")

%% Hluti 3 
%dæmi 1
clc
y0=[0.8; 0.1; 8];
n=1000;T=400; b1=3;
t=linspace(0,T,n+1);
rk=kutta3(y0,T,n,b1);
plot(t, rk(1, :), 'r', t, rk(2, :), 'b', t, rk(3, :), "g", 'LineWidth', 1);
legend("X plöntur", "Y plöntuætur", "Z kjötætur")

% Dæmi 2
figure
plot3(rk(1, :),rk(2, :), rk(3,:))

%% Dæmi 2 aftur

%% Dæmi 2 aftur
y0_1=[0.8; 0.1; 8];
y0_2=[1;0.5;5];
y0_3=[4;2;4];
n=100000;T=1000; b1=3;
t=linspace(0,T,n+1);
rk_1=kutta3(y0_1,T,n,b1);
rk_2=kutta3(y0_2,T,n,b1);
rk_3=kutta3(y0_3,T,n,b1);
t=linspace(0,T,n+1);
hold on

plot(t, rk_1(1, :), 'b', t, rk_1(2, :),"--b", t , rk_1(3, :), '--b', 'LineWidth', 1.5);
plot(t, rk_2(1, :), 'g', t, rk_2(2, :), '--g',t , rk_2(3, :), '--g', 'LineWidth', 1.5);
plot(t, rk_3(1, :), 'r', t, rk_3(2, :), '--r',t , rk_3(3, :), '--r', 'LineWidth', 1.5);
legend("y0_1","y0_2","y0_3")
title("Mismunandi upphafsgildi")
%Greinilegt að Allt stefnir á sama stað þegar T er nógu stórt. 
%% dæmi 3

b1=[0.5, 3];
y0=[0.8; 0.1; 8];
b1=linspace(0.5,3,6);
n=1000;T=1000;
t=linspace(0,T,n+1);
for i=1:length(b1)
    rk=kutta3(y0,T,n,b1(i));
    figure
    plot(t, rk(1, :), 'r', t, rk(2, :), 'b', t, rk(3, :), "g", 'LineWidth', 1);
    title("b1= ", b1(i))
    legend("X plöntur", "Y plöntuætur", "Z kjötætur")
end


%% Dæmi 3 aftur

y0=[0.8; 0.1; 8];
b1=linspace(0.5,3,6);
n=10000;T=1000;
t=linspace(0,T,n+1);
figure
names=[];
for i=1:length(b1)
    rk=kutta3(y0,T,n,b1(i));
    hold on
    names= [ names , "z kjötætur, b=" + b1(i)];
    plot(t, rk(3, :), 'LineWidth', 1,'DisplayName', names(i))
    hold off
end
title("Z (Kjötætur) miðað við mimsunandi b_1 fasta")
legend("show")

figure
names=[];
b1=[0.5 1 2 3];
T=1000;
t=linspace(0,T,n+1);
for i=1:length(b1)
    rk=kutta3(y0,T,n,b1(i));
    hold on
    plot(t, rk(1, :),"g", t, rk(2, :), "r",'LineWidth', 1)
    hold off
end
title("plöntur og plöntuætur fyrir mismunandi b")
legend("show")


