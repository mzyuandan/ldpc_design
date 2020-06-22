clc;
clear all;
close all;

mb = 12;
nb = 18;
dv = 7;
dc = 5;
evl_num = 50;
BER = 10^(-5);
Qval_t = 2*(qfuncinv(BER)^2);
test_num = 1000;

threshold = zeros(1, test_num);
SNR_t = zeros(1, test_num);
var_vect = zeros(test_num, dv-1);
chk_vect = zeros(test_num, 3);
tt = 1;
while tt <= test_num

err = 1;
while err==1
[var_vect, edg_num] = gen_var_vect(mb, nb, dv);
[chk_vect, err] = gen_chk_vect(mb, edg_num, dc);
end

lam = zeros(1, dv-1);
for jj=2:1:dv
    lam(jj-1) = jj*var_vect(jj-1) / edg_num;
end

rou = zeros(1, dc-1);
for ii=2:1:dc-3
    rou(ii-1) = 0;
end
for ii=dc-2:1:dc
    rou(ii-1) = ii*chk_vect(ii-dc+3) / edg_num;
end

threshold = threshold_cal_ir_new( lam, rou, evl_num, Qval_t );

SNR_t = 10*log10(1/(2*(threshold^2)));

if (SNR_t<=-4)
    var_vect_good(tt,:) = var_vect;
    chk_vect_good(tt,:) = chk_vect;
    thold(tt) = threshold;
    SNR(tt) = SNR_t;
    fprintf('tt=%d, SNR=%f\n', tt, SNR(tt));
    tt = tt+1;
end
end
