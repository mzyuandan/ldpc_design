function [ threshold ] = threshold_cal( dc, dv, evl_num )
%UNTITLED1 Summary of this function goes here
%   Detailed explanation goes here


% clc
% clear all
% dc = 3;
% dv =4;
% evl_num = 1000;

threshold = 0.0;

for ii=1:1:10000
    threshold = threshold + 0.001;
    [mu, mv] = evl(threshold, dc, dv, evl_num);
    if (mu(evl_num) < 10.634300)
        break
    end
end
end

