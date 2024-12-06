% Hreyfi.m

clear all
close all
L=30;lambda=0.8;EI=1.09*10^(10);
k=0.1565;
omega=2.8578*10^3;
delta=-3.5;

T = 2.0*pi; % Ein lota

% Búa til glugga til að teikna hreyfimyndinna í
% Við felum hann á meðan við teiknum.
hfig = figure('visible', 'off');

FPS = 30; % Fjöldi ramma á sek í hreyfimyndinni.
          % Til að hægja á hreyfimyndinni er hægt að auka FPS
          % og halda FPS_PLAY föstu.

FPS_PLAY = 30; % Fjöldi ramma á sek þegar við spilum hreyfimyndinna.
               % Best að hafa það 30.

t_start = 0.0; % Byrjunar tími í sek
t_end = 2*T; % Enda tími í sek (Tvær lotur)

% Reikna út fjölda ramma sem þarf og stærð eins tíma skrefs.
n = ceil(FPS*(t_end - t_start));
dt = (t_end - t_start)/n;

F(1:n) = getframe(hfig);
t = zeros(1, n);

% Teikna hvert tíma skref
for i = 1:n
    % Hreinsa allt úr glugganum
    cla
    hold on
    
    % Stilla stærð ása.
    % Best er að þeir séu alltaf jafn stórir í hverjum ramma.
    axis([0.0, L, -10.0, 10.0]) %[x_min, x_max, y_min, y_max]
    %axis equal
    
    % Reikna tímann
    t(i) = t_start + i*dt;
    
    % Teikna dót - þessu þarf að breyta og finna (x,y) hnit
    x=linspace(0,L,100);
    plot(x, A(x,k,L)*B(t(i),k,delta,L,omega))
    
    % Ef við viljum sýna tímann á grafinu
    text(5.0, 5.0, strcat(num2str(t(i)), ' s'))
    
    grid on % Grid on eða off
    
    % Merkja ásanna og setja titil
    xlabel('x [eining]')
    ylabel('y [eining]')
    title('Hreyfimynd af sin')
    
    % Geyma rammann
    F(i) = getframe(hfig);
end

% Sýna hreyfimyndina einu sinni
close all
movie(gcf, F, 1, FPS_PLAY)

% Skrifa hreyfimyndinni niður a diskinn
v = VideoWriter('Video_1.mp4', 'MPEG-4');
v.FrameRate = FPS_PLAY;
open(v);
writeVideo(v, F);
close(v);