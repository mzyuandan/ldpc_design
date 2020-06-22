function [ mu, mv ] = evl( sigma, dv, dc, evl_num )
%DESITY_EVL_G Summary of this function goes here
%   Detailed explanation goes here

% clc
% claer all
% sigma=1.252;
% dv=3;
% dc=4;
% evl_num=1000;

mu = zeros(1, evl_num);
mv = zeros(1, evl_num);
m0 = 2/(sigma*sigma);
mu(1) = inverse_cal(1-((1-base_cal(m0))^(dc-1)));
for l=2:1:evl_num  
     mu(l) = inverse_cal(1-((1-base_cal(m0+((dv-1)*mu(l-1))))^(dc-1)));
     mv(l) = m0+(dv-1)*mu(l);
end
end

