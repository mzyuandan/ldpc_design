function [ u ] = inverse_cal( v )
%INVERSE_CAL Summary of this function goes here
%   Detailed explanation goes here

load u_tab.mat;
if (v==0)
    u = 1000;
elseif (v>=0.0385 && v<=1.022)
    u = ((0.0218-log(v))/0.4517)^(1/0.86);
elseif (v<0.0385 && log10(v)>=-110)
    u = u_tab(round(-1000*(log10(v)+1.404)));
elseif (log10(v)>-110)
    u = 1000;
end

end
