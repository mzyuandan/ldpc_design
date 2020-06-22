function [ threshold ] = threshold_cal_ir_new( lam, rou, evl_num, Qval_t, th_step, th_start )
%DIST_TEST Summary of this function goes here
%   Detailed explanation goes here

% clc
% clear all
% dc = 3;
% dv =4;
% evl_num = 1000;

threshold = th_start;

for ii=1:1:10000
    threshold = threshold + th_step*10;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end
threshold = threshold - th_step;
for ii=1:1:10000
    threshold = threshold + th_step*10;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end
threshold = threshold - th_step;
for ii=1:1:10000
    threshold = threshold + th_step;
    [mv_all, mu, mv, m0] = evl_ir_new(threshold, lam, rou, evl_num);
    if (mv_all < Qval_t)
        break
    end
end

end


