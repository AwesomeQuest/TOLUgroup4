% Hreyfi.m

clear all
close all
L=30;lambda=0.8;EI=1.09*10^(10);
k=0.1565;
omega=2.8578*10^3;
delta=-3.5;

T = 2.0*pi; % Ein lota

% B�a til glugga til a� teikna hreyfimyndinna �
% Vi� felum hann � me�an vi� teiknum.
hfig = figure('visible', 'off');

FPS = 30; % Fj�ldi ramma � sek � hreyfimyndinni.
          % Til a� h�gja � hreyfimyndinni er h�gt a� auka FPS
          % og halda FPS_PLAY f�stu.

FPS_PLAY = 30; % Fj�ldi ramma � sek �egar vi� spilum hreyfimyndinna.
               % Best a� hafa �a� 30.

t_start = 0.0; % Byrjunar t�mi � sek
t_end = 2*T; % Enda t�mi � sek (Tv�r lotur)

% Reikna �t fj�lda ramma sem �arf og st�r� eins t�ma skrefs.
n = ceil(FPS*(t_end - t_start));
dt = (t_end - t_start)/n;

F(1:n) = getframe(hfig);
t = zeros(1, n);

% Teikna hvert t�ma skref
for i = 1:n
    % Hreinsa allt �r glugganum
    cla
    hold on
    
    % Stilla st�r� �sa.
    % Best er a� �eir s�u alltaf jafn st�rir � hverjum ramma.
    axis([0.0, L, -10.0, 10.0]) %[x_min, x_max, y_min, y_max]
    %axis equal
    
    % Reikna t�mann
    t(i) = t_start + i*dt;
    
    % Teikna d�t - �essu �arf a� breyta og finna (x,y) hnit
    x=linspace(0,L,100);
    plot(x, A(x,k,L)*B(t(i),k,delta,L,omega))
    
    % Ef vi� viljum s�na t�mann � grafinu
    text(5.0, 5.0, strcat(num2str(t(i)), ' s'))
    
    grid on % Grid on e�a off
    
    % Merkja �sanna og setja titil
    xlabel('x [eining]')
    ylabel('y [eining]')
    title('Hreyfimynd af sin')
    
    % Geyma rammann
    F(i) = getframe(hfig);
end

% S�na hreyfimyndina einu sinni
close all
movie(gcf, F, 1, FPS_PLAY)

% Skrifa hreyfimyndinni ni�ur a diskinn
v = VideoWriter('Video_1.mp4', 'MPEG-4');
v.FrameRate = FPS_PLAY;
open(v);
writeVideo(v, F);
close(v);