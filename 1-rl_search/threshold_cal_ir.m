function [ threshold ] = threshold_cal_ir( lam, rou, evl_num )
%THRESHOLD_CAL_IR Summary of this function goes here
%   Detailed explanation goes here

% clc
% clear all
% dc = 3;
% dv =4;
% evl_num = 1000;

threshold = 0.0;

for ii=1:1:10000
    threshold = threshold + 0.001;
    [mu, mv] = evl_ir(threshold, lam, rou, evl_num);
    if (mu(evl_num) < 10.0)
        break
    end
end

end
