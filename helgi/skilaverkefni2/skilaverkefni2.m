%% Daemi 1
clc; clear; close all;

% Tokum gognin ur toflunni sem gefin er i verkefnalysingunni
T = (2:10) .*10 + 273.15; % Breyta Celsius i Kelvin
T = T';
mu_data = [1.003, 0.799, 0.657, 0.548, 0.467, 0.405, 0.355, 0.316, 0.283]'*1e-3;

% Gogn til thess ad plotta!
T_plot = linspace(20,100);

% Skilgreinum dalka fylkisins A
A1 = ones(size(mu_data));
A2 = T.^(-1);
A3 = T.^(-2);

% Setjum i fylkid A
A = [A1,A2,A3];

% Leysum jofnuhneppid A*beta = mu
beta_1 = A\log(mu_data);

% Reiknum Root Mean Square Error (RMSE) af nalguninni okkar
error_1 = norm(A*beta_1 - log(mu_data));

% set(gca,'TickLabelInterpreter','latex');
% Plottum graf af seigd sem fall af hitastigi (Celcius fyrir einfaldleika)
plot(T - 273.15,mu_data,".",MarkerSize=15);

% Plottum einnig linulegu nalgunina a sama graf
% Notum fallid mu_func.m sem eg bjo til ut fra fallinu fyrir seigd i verkefnalysingunni
hold on
plot(T - 273.15,mu_func1(T,beta_1),LineWidth=1.5);
xlabel("Hitastig [\circC]");
ylabel("Seigð [10^{-3} Pa\cdots]")
grid on
title("Dæmi 1 - Aðferð minnstu ferninga")

%% Daemi 2

% Nota föllin:
% - gaussnewton.m
% - F.m
% - jac.m

x0 = zeros(3,1);
tol = 1e-5;
beta_2 = gaussnewton(x0,tol);

% Metum hvor gefur meiri skekkju
error_2 = norm(F(beta_2));
error_1_re = norm(F(beta_1));

if error_1 < error_2
    disp("Adferd minnstu ferninga skilar minni skekkju!")
elseif error_1 == error_2
    disp("Badar adferdir skila somu skekkju!")
else
    disp("Adferd Gauss-Newton skilar minni skekkju!")
end

% set(gca,'TickLabelInterpreter','latex');
% Plottum graf af seigd sem fall af hitastigi (Celcius fyrir einfaldleika)
plot(T - 273.15,mu_data,".",MarkerSize=15);

% Plottum einnig linulegu nalgunina a sama graf
% Notum fallid mu_func.m sem eg bjo til ut fra fallinu fyrir seigd í verkefnalysingunni
hold on
plot(T - 273.15,mu_func1(T,beta_2),LineWidth=1.5);
xlabel("Hitastig [\circC]");
ylabel("Seigð [10^{-3} Pa\cdots]");
grid on
title("Dæmi 2 - Aðferð Gauss-Newton");

%% Daemi 3

% Adferd minnstu ferninga
% Skilgreinum dalka fylkisins B
B1 = T.^(-1);
B2 = T;
B3 = T.^(2);

% Setjum i fylkid B
B = [B1,B2,B3];

% Leysum jofnuhneppid A*beta = mu
beta_3 = B\log(mu_data);

% Reiknum Root Mean Square Error af nalguninni okkar
error_3 = norm(B*beta_3 - log(mu_data));

% set(gca,'TickLabelInterpreter','latex');
% Plottum graf af seigd sem fall af hitastigi (Celcius fyrir einfaldleika)
plot(T - 273.15,mu_data,".",MarkerSize=15);

% Plottum einnig linulegu nalgunina a sama graf
% Notum fallid mu_func.m sem eg bjo til ut fra fallinu fyrir seigd í verkefnalysingunni
hold on
plot(T - 273.15,mu_func2(T,beta_3),LineWidth=1.5);
xlabel("Hitastig [\circC]");
ylabel("Seigð [10^{-3} Pa\cdots]");
grid on
title("Dæmi 3 - Aðferð minnstu ferninga");

% Adferd Gauss-Newton
beta_4 = gaussnewton2(x0,tol);

figure();
plot(T - 273.15,mu_data,".",MarkerSize=15);
hold on
plot(T - 273.15,mu_func2(T,beta_4),LineWidth=1.5);
xlabel("Hitastig [\circC]");
ylabel("Seigð [10^{-3} Pa\cdots]");
grid on
title("Dæmi 3 - Aðferð Gauss Newton");

% Metum hvor gefur meiri skekkju
error_4 = norm(F2(beta_4));
error_3_re = norm(F2(beta_3));

%% Daemi 4

n = 2; 
L = 0.1; % lengd milli platna i m
V0 = 0.01; % hradi efri plotu i m/s

% Notum nalgun a mu fallinu sem gefur okkur minnstu skekkjuna
% Thad er Gauss-Newton ur daemi 2

K = V0 / simpson(0,L,n);

% Buum til vigur med nokkrum y gildum fra 0 uppi L og annan fyrir hradann (gogn)
y = linspace(0,L,100);
V_data1 = zeros(size(y));

% Lykkjan her setur hvert stak vigursins y inn i jofnuna V og geymir
% hradann i hradavigrinum
for i = 1:length(y)
    V_data1(i) = K * simpson(0,y(i),n);
end

% Plottum gognin a milli 0 og L
plot(y,V_data1,LineWidth=1.5);
xlabel("Bil [m]");
ylabel("Hraði vökva [m/s]");
title("V(y)");

%% Daemi 5

% Skilgreinum vigur med mismunandi n-um og tomann gagnavigur
n_vector = [2 4 8 16 32 64 128 256];
V_data2 = zeros(size(n_vector));
% Skodum y = L/2
dist1 = L/2;

% Veljum n = 2^14 til thess ad mida vid (mjooooog nakvaemt)
V_real = V0 * (simpson(0,dist1,n^14) / simpson(0,L,n^14));

% Setjum inn hvert n fyrir sama y og sofnum i gagnavigurinn
for i = 1:length(n_vector)
    V_data2(i) = V0 * simpson(0,dist1,n_vector(i)) / simpson(0,L,n_vector(i));
end

V_error = abs(V_data2 - V_real);

log2(abs(diff(V_data2)));

plot(log2(n_vector),V_error,".",MarkerSize=15);
xt = get(gca, 'XTick');
set(gca, 'XTickLabel', 2.^xt);

%% Daemi 6

tol_adap = 0.5e-10;

% Trapisuadferd
K1 = V0 / adapquad_trap(0,L,tol_adap);
V1 = K1 * adapquad_trap(0,dist1,tol_adap);

% Adferd Simpson
K2 = V0 / adapquad(0,L,tol_adap);
V2 = K2 * adapquad(0,dist1,tol_adap);

%% Daemi 7

% Bua til glugga til ad teikna hreyfimyndinna i
% Vid felum hann a medan vid teiknum.
hfig = figure('visible', 'off');

FPS = 30; % Fjoldi ramma a sek i hreyfimyndinni.
          % Til ad haegja a hreyfimyndinni er haegt ad auka FPS
          % og halda FPS_PLAY fostu.

FPS_PLAY = 30; % Fjoldi ramma a sek thegar vid spilum hreyfimyndinna.
               % Best ad hafa thad 30.

t_start = 0.0; % Byrjunartimi i sek
t_end = 10; 

% Reikna ut fjolda ramma sem tharf og staerd eins timaskrefs.
n = ceil(FPS*(t_end - t_start));
dt = (t_end - t_start)/n;

Frames(1:n) = getframe(hfig);
t = zeros(1, n);

% Teikna hvert timaskref
for i = 1:n
    % Hreinsa allt ur glugganum
    cla
    hold on
    
    % Stilla staerd asa.
    % Best er ad their seu alltaf jafn storir i hverjum ramma.
    axis([0.0, L, 0, V0]) %[x_min, x_max, y_min, y_max]
    %axis equal
    
    % Reikna timann
    t(i) = t_start + i*dt;
    
    % Teikna dot - Thessu tharf ad breyta og finna (x,y) hnit
%     x=linspace(0,L,100);
%     plot(x, A(x,k,L)*B(t(i),k,delta,L,omega))
    plot(y, V_data1)
    
    % Ef vid viljum syna timann a grafinu
    text(5.0, 5.0, strcat(num2str(t(i)), ' s'))
    
    grid on % Grid on eda off
    
    % Merkja asanna og setja titil
    xlabel('Vegalengd [m]')
    ylabel('Hraði [m/s]')
    title(['Hreyfimynd af lóðréttri línu í vökvanum'])
    
    % Geyma rammann
    Frames(i) = getframe(hfig);
end

% Syna hreyfimyndina einu sinni
close all
movie(gcf, Frames, 1, FPS_PLAY)

% Skrifa hreyfimyndinni nidur a diskinn
v = VideoWriter('Video_2.mp4', 'MPEG-4');
v.FrameRate = FPS_PLAY;
open(v);
writeVideo(v, Frames);
close(v);

%% Daemi 8

% Daemi 4 endurtekid
% Adeins tharf ad breyta mu.m skranni
V_data3 = zeros(size(y));

% Lykkjan her setur hvert stak vigursins y inn i jofnuna V og geymir
% hradann i hradavigrinum
for i = 1:length(y)
    V_data3(i) = K * simpson(0,y(i),n);
end
