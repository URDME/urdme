
%Chrishani Jayaweera 2020-11-04
%Euler forward 

%Define start values
%Define step size

function [t, y] = Euler(t0, y0, h, tn)
    t = (t0:h:tn)';
    y = zeros(size(t));
    y(1) = y0;
    for i = 1:1:length(t) - 1
        y(i + 1) = y(i) + h * f(y(i), t(i));
        
        
        