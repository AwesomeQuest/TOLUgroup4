clc; clear; close all;
y0 = [0.8; 0.1; 8];
b1 = 0.5:0.5:3;
n = 40000; T = 900;
figure;
names = strings(1,length(b1));
x = 100;
y = 200;
height = 400;
width = 1000;

for i = 1:length(b1)
    [time, yvec_rk_z] = RKsolver3(y0,n,T,b1(i));
    hold on
    names(i) = "Z Kjötætur, b_1=" ;
    plot(time, yvec_rk_z(3, :), 'LineWidth', 1,'DisplayName', names(i)+ b1(i))
    hold off
end

title("Z Kjötætur fyrir mismunandi b_1 fasta")
xlabel("Tími");
ylabel('Kjötætur');
legend("show",'location', 'northwest')

filename = sprintf('h3_d3_Z.svg');
% saveas(gcf, filename);

figure;
T = 1000;
t = linspace(0,T,n+1);
names = strings(1,length(b1));

for i = 1:length(b1)
    [timev, yvec_rk_x] = RKsolver3(y0,n,T,b1(i));
    hold on
    names(i) = "X Plöntur, b_1=" ;
    plot(timev, yvec_rk_x(1, :),'LineWidth', 1,'DisplayName', names(i)+ b1(i))
end

title("X Plöntur fyrir mismunandi b_1 fasta")
xlabel("Tími");
ylabel('Plöntur');
legend('show')
set(gcf,'position',[x,y,width,height])

filename = sprintf('h3_d3_X.svg');
% saveas(gcf, filename);

figure;
names = strings(1,length(b1));
for i = 1:length(b1)
    [timev, yvec_rk_y] = RKsolver3(y0,n,T,b1(i));
    hold on
    names(i) = "Y Plöntuætur, b_1=" ;
    plot(timev, yvec_rk_y(2, :), 'lineWidth',1,'DisplayName', names(i)+ b1(i))
end

hold off
title("Y Plöntuætur fyrir mismunandi b_1")
xlabel("Tími");
ylabel('Plöntuætur');
legend("show")
set(gcf,'position',[x,y,width,height])

filename = sprintf('h3_d3_Y.svg');
% saveas(gcf, filename);

names = strings(1,length(b1));
for i = 1:length(b1)
    [timev, yvec_rk_xyz]=RKsolver3(y0,n,T,b1(i));
    names(i) = "b_1=" ;
    figure;
    plot3(yvec_rk_xyz(1,:),yvec_rk_xyz(2,:),yvec_rk_xyz(3,:),"LineWidth",0.8,'DisplayName', names(i)+ b1(i))
    legend('show','location', 'west')
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    filename = sprintf('h3_d3_%.2f.svg', b1(i));
%     saveas(gcf, filename);
end

hold off 