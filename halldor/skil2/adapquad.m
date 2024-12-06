function x=adapquad(a,b,tol, interval_type)

persistent nefnari_intervals teljari_intervals
if isempty(nefnari_intervals)
    nefnari_intervals = zeros(1, 2); 
end
if isempty(teljari_intervals)
    teljari_intervals = zeros(1, 2); 
end

% displayar hversu mörg intervals eru
if nargin == 0
    fprintf('Numerator intervals: %d\n', size(nefnari_intervals, 1));
    fprintf('Denominator intervals: %d\n', size(teljari_intervals, 1));
    fprintf('Total intervals: %d\n', size(nefnari_intervals, 1) + size(teljari_intervals, 1));
    nefnari_intervals = zeros(1, 2); 
    teljari_intervals = zeros(1, 2); 
    x = [];
    return
end

c=(a+b)/2;
sab=simpson(a,b,1);sac=simpson(a,c,1);scb=simpson(c,b,1);

% athugar fyrir duplicates
if strcmp(interval_type, 'numerator')
    if ~ismember([a, b], nefnari_intervals, 'rows')
        nefnari_intervals = [nefnari_intervals; a, b]; 
    end
elseif strcmp(interval_type, 'denominator')
    if ~ismember([a, b], teljari_intervals, 'rows')
        teljari_intervals = [teljari_intervals; a, b];
    end
end

if abs(sab-sac-scb)<10*tol
    x=sac+scb;
else
    x1 = adapquad(a, c, tol / 2, interval_type);
    x2 = adapquad(c, b, tol / 2, interval_type);
    x = x1 + x2;
end
end