%% Verkefni 4 allt saman
%dæmi 4 TEKUR ~10 SEK að keyra
%tekið úr comporgbetter.m
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

%% Dæmi 4 myndband tekur 20 sek að Keyra kóðan ----------------------------------------------------
clc
T=30;N=100;M=N;L=5;
D = 0.01;v = 0.1;h = L / M;k = T / N;
[w,ts,xs] = iterdiffv2(T,N,M);
max_conc = max(w(:));
FPS = 30;                 % Frames per second for playback
n_frames = N + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames
x = linspace(0, L, M+1);
for j = 1:n_frames
    % Clear the figure
    clf;
    hold on;
    
    % Plot pollutant concentration
    plot(x, w(:, j), 'b-', 'LineWidth', 2);
    
    % Set axis limits and labels
    axis([0, L, 0, max_conc * 1.1]); % Adjust limits based on expected results
    xlabel('Position (m)');
    ylabel('Concentration (kg/m^3)');
    title(sprintf('Pollutant Distribution at t = %.2f s', (j-1) * k));
    grid on;
    
    % Store frame
    F(j) = getframe(hfig);
end

% Play animation
close all;
movie(gcf, F, 1, FPS);

% % Save animation as video
% v = VideoWriter('Pollutant_Animation2.mp4', 'MPEG-4');
% v.FrameRate = FPS;
% open(v);
% writeVideo(v, F);
% close(v);
% 
% [~, max_idx] = max(w(:, end));
% disp(['Maximum concentration at t = 30 s: ', num2str(max(w(:, end))), ' at x = ', num2str(x(max_idx))]);



%% Dæmi 5 Tekur 5 sekúndur
% Parameters for Part 5
close all
T=30;N1=100;M1=N1;
N2=1000;M2=N2;
k1 = T / N1;
k2= T / N2;
[w1,t1,x1] = iterdiffv1(T,N1,M1);
[w2,t2,x2] = iterdiffv1(T,N2,M2);
time_steps = [0, T/6, 2*T/6, 3*T/6, 4*T/6, 5*T/6];
% Plot results for comparison
for t_idx = 1:length(time_steps)
    t_plot = time_steps(t_idx);
    idx1 = round(t_plot / k1) + 1; % Find index for Part 4
    idx2 = round(t_plot / k2) + 1; % Find index for Part 5
    
    figure;
    hold on;
    plot(x1, w1(:, idx1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dæmi 4 (m=n=100)');
    plot(x2, w2(:, idx2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Dæmi 5 (m=n=1000)');
    hold off;
    
    xlabel('Staðsetning (m)');
    ylabel('Þéttleiki (kg/m^3)');
    title(sprintf('Dreifing mengunarefninsin við t = %.0f/6T', t_idx-1));
    legend('Location', 'best');
    grid on;
end

%% Dæmi 6 ---------------------------------------------------------------
%% Viðmiðunargildi f dæmi 6 
% tekur 10 sekúndur að keyra
n_real=10000;T=30;L=5;
[mat_real,~,~] = iterdiffv2(T,n_real,n_real);
%% Dæmi 6 t=30 loglog plots  mynd 4
%tekur 15 sekúndur að keyra
close all; clc
T=30;L=5;
n_high=3000;
skekkja=zeros(2,n_high/100);
n=linspace(100,n_high,n_high/100);
n=round(n,0);
real_high_val=max(mat_real(:,end));
count=1;
for i=n
    [mat1,ts,x_1] = iterdiffv1(T,i,i); %fall sem notar upprunalega A fylki og b
    [mat2,td,x_2] = iterdiffv2(T,i,i); %fall sem skiptir um aðferð fyrir t eftir fyrstu ítrun

    high_val1=max(mat1(:,end)); %finnur hæsta gildið
    high_val2=max(mat2(:,end));

    skekkja1=abs(real_high_val-high_val1)/(real_high_val)*100; %reiknar skekkju 
    skekkja2=abs(real_high_val-high_val2)/(real_high_val)*100;

    skekkja(1,count)=skekkja1;
    skekkja(2,count)=skekkja2;
    count=count+1;
end

% plot(n,skekkja)
% legend("show")
% plot(x_1,mat1(:,end))
% hold on
% plot(x_2,mat2(:,end))
% hold off
%% mynd 4 dæmi 6
% tekur 1 sek að keyra
logn=log(n);log_skekkja=log(skekkja);
p1 = polyfit(logn,log_skekkja(1,:),1);
p2 = polyfit(logn,log_skekkja(2,:),1);
% disp(p1(1))
% disp(p2(1))
figure;
loglog(n,skekkja(1,:),'r',LineWidth=2)
txt1="Hallatala = "+ p1(1) ;text(600,2,txt1,"FontSize",14) 
xlabel("log(n)","FontSize",14)
ylabel("log(skekkja)","FontSize",14)
% title("loglog plot af sekkkju og n")
saveas(gcf,"d6_loglog1.svg")

figure;
loglog(n,skekkja(2,:),"b",LineWidth=2)
txt1="Hallatala = "+ p2(1) ;text(600,0.02,txt1,'FontSize',14) 
xlabel("log(n)","FontSize",14)
ylabel("log(skekkja)","FontSize",14)
% title("loglog plot af sekkkju og n")
saveas(gcf,"d6_loglog2.svg")


%% Dæmi 6
% tekur 10 sekúndur að keyra
n_real=10000;T=30;
[mat_real,~,~] = iterdiffv1(T,n_real,n_real);
%% mynd 5.a
%tekur 80 sekúndur að keyra
T=30;
n_high=9000;
skekkja=zeros(1,n_high/300);
n=linspace(300,n_high,n_high/300);
n=round(n,0);
real_high_val=max(mat_real(:,end));
count=1;
for i=n
    [mat1,ts,x_1] = iterdiffv1(T,i,i); %fall sem notar upprunalega A fylki og b

    high_val1=max(mat1(:,end)); %finnur hæsta gildið

    skekkja1=abs(real_high_val-high_val1)/(real_high_val)*100; %reiknar skekkju 

    skekkja(1,count)=skekkja1;
    count=count+1;
end
%%  mynd 5.a
% Tekur 1 sekúndu að keyra
plot(n,skekkja,'r',LineWidth=1)
set(gca,'fontsize',14)
yline(0.1,'g', LineWidth=1)
xline(4200,'b',LineWidth=1)
xlabel("n")
ylabel("Skekkja [%]")

%% Raungildi 2.stigs skekkja metin á 1 stigs
% mynd 5.b
% Tekur 80-90 skúndur að keyra
n_real=100;T=30;
[mat_real,~,~] = iterdiffv2(T,n_real,n_real);

T=30;L=5;
n_high=9000;
skekkja=zeros(1,n_high/300);
n=linspace(300,n_high,n_high/300);
n=round(n,0);
real_high_val=max(mat_real(:,end));
count=1;
for i=n
    [mat1,ts,x_1] = iterdiffv1(T,i,i); %fall sem notar upprunalega A fylki og b

    high_val1=max(mat1(:,end)); %finnur hæsta gildið

    skekkja1=abs(real_high_val-high_val1)/(real_high_val)*100; %reiknar skekkju 

    skekkja(1,count)=skekkja1;
    count=count+1;
end
%% Hér er verið að kanna skekkjuna þear n=4000  og raungildi er reiknað m.v. n=100 en 2. stigs aðferð
n=40000;
[mat1,ts,x_1] = iterdiffv1(T,n,n); %fall sem notar upprunalega A fylki og b
high_val1=max(mat1(:,end)); %finnur hæsta gildið
skekkja1=abs(real_high_val-high_val1)/(real_high_val)*100 %reiknar skekkju 


%% mynd 5.b
% tekur 1 sek 

plot(n,skekkja,'r',LineWidth=1)
yline(0.1,'g', LineWidth=1)
% xline(2600,'b',LineWidth=1)
set(gca,'fontsize',14)
xlabel("n")
ylabel("Skekkja [%]")

%% Raungildi 2 stigs og skekkja metin á 2 stigs
% mynd 5.c
%Tekur 3 sekúndur að keyra
clear,clc
n_real=4000;T=30;
[mat_real,~,~] = iterdiffv2(T,n_real,n_real);
T=30;L=5;
n_high=500;
% n=linspace(100,n_high,n_high/10);
n=100:10:500;
skekkja=zeros(1,length(n));
real_high_val=max(mat_real(:,end));
count=1;
for i=n
    [mat1,ts,x_1] = iterdiffv2(T,i,i); %fall sem notar upprunalega A fylki og b
    high_val1=max(mat1(:,end)); %finnur hæsta gildið

    skekkja1=abs(real_high_val-high_val1)/(real_high_val)*100; %reiknar skekkju 

    skekkja(1,count)=skekkja1;
    count=count+1;
end
%% mynd 5.c
plot(n,skekkja,'r',LineWidth=1)
yline(0.1,'g', LineWidth=1)
% xline(2600,'b',LineWidth=1)
xlabel("n")
ylabel("Skekkja [%]")
set(gca,'fontsize',14)
xline(120)

%% Dæmi 8----------------------------------------------------
clc;clear;close all;
T = 100;
N = 100;
M = N;
k = T/N;


[us1,ts1,xs1] = iterdiffd8(T,N,M);
% [us1,ts,xs1] = iterdiffd8org(T,N,M);
% plot(xs1,us1(:,1))

for i = 1:6
    figure;
    j = -i + 7
    plot(xs1,us1(:,round(N/6*j)+1),LineWidth=2)
    ylim([0 1])
    xlabel("Staðsetning [m]")
    ylabel("Þéttleiki [kg/m^3]")
    grid on
    txt = "d8plot" + j + ".svg"
    saveas(gcf,txt)
end


% plot3(xs1,ts1,us1);
%% LETS ANIMATE THIS MF

% Animation parameters
FPS = 15;                 % Frames per second for playback
n_frames = N + 1;         % Total number of frames
hfig = figure('visible', 'off'); % Hidden figure window for rendering
F(1:n_frames) = getframe(hfig); % Preallocate frames

for i = 1:N+1
    clf;
    hold on
    plot(xs1,us1(:,i),LineWidth=2)
%     hold on
%     plot(xs2,us2(:,i))
%     plot(xs3,us3(:,10*i))
%     legend("original", "modified", "many iterations original")
    ylim([0,1])
%     pause(0.1)
    xlabel('Staðsetning [m]');
    ylabel('Þéttleiki [kg/m^3]');
    title(sprintf('Dreifing mengunarefnisins við t = %.2f s', (i-1) * k));
    grid on
    F(i) = getframe(hfig);
%     hold off
end

% Play animation
close all;
movie(gcf, F, 1, FPS);

% Save animation as video
v = VideoWriter('Pollutant_Animation_d8.mp4', 'MPEG-4');
v.FrameRate = FPS;
open(v);
writeVideo(v, F);
close(v);

%% Dæmi 9 -----------------------------------------------------------------
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
