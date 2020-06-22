function [ mu, mv, m0 ] = evl_ir( sigma, lam, rou, evl_num )
%EVL_IR Summary of this function goes here
%   Detailed explanation goes here

% clc
% clear all
% sigma=1.251;
% lam=[0 0 1];
% rou=[0 0 0 1];
% evl_num=1000;

[tmp,dl] = size(lam);
[tmp,dr] = size(rou);
mu = zeros(1, evl_num);
mv = zeros(1, evl_num);
m0 = 2/(sigma*sigma);
v2u = 0;
for jj=2:1:dl
    v2u = v2u + lam(jj)*base_cal(m0);
end
for ii=2:1:dr
    mu(1) = mu(1) + rou(ii)*inverse_cal(1-(1-v2u)^(ii-1));
end
for jj=2:1:dl
    mv(1) = mv(1) + lam(jj)*(m0+(jj-1)*mu(1));
end
    
for l=2:1:evl_num
    v2u = 0.0;
    for jj=2:1:dl
        v2u = v2u + lam(jj)*base_cal(m0+(jj-1)*mu(l-1));
    end
    for ii=2:1:dr
        mu(l) = mu(l) + rou(ii)*inverse_cal(1-(1-v2u)^(ii-1));
    end
    for jj=2:1:dl
    mv(l) = mv(l) + lam(jj)*(m0+(jj-1)*mu(l));
end
end
%end
