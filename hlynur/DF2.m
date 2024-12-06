function y=DF2(y,w)
    y1=y(1);
    y2=y(2);
    r1 = 1; r2 = 0.1; k = 7; d = 1; j = 1;;
    dy1=r1*y1*(1-y1/k)-w*y1*y2/(d+y1); 
    dy2=r2*y2*(1-(j*y2)/y1);
    y=[dy1;dy2];
end
