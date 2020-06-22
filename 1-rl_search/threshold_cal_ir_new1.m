function [ threshold ] = threshold_cal_ir_new1( lam, rou, evl_num, Qval_t, th_step, th_start )
%DIST_TEST Summary of this function goes here
%   Detailed explanation goes here

% clc
% clear all
% dc = 3;
% dv =4;
% evl_num = 1000;

threshold = th_start;
th_stepX100 = th_step * 100;
th_stepX10 = th_step * 10;

for ii=1:1:10000
    threshold = threshold + th_stepX100;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end
threshold = threshold - th_stepX100;
for ii=1:1:10000
    threshold = threshold + th_stepX10;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end
threshold = threshold - th_stepX10;
for ii=1:1:10000
    threshold = threshold + th_step;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end
threshold = threshold - th_step;
end


