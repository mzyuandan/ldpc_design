%function [ vmax_mu, vmax_mv ] = quantiz_bit( lam, rou, evl_num, arcq, code_rate )                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  )
%QUANTIZ_BIT Summary of this function goes here
%   Detailed explanation goes here

clc
% clear all
% lam=[0 0 1];
% rou=[0 0 0 1];
evl_num=20;
grade = 0;
arcq = 0.0;
code_rate = 0.5;
satu = 11;
 
snr_start = 0.0;
snr_step = 0.2;
snr_end = 2.2;

fp = fopen('vmax_c.txt', 'w');
fp1 = fopen('vmax_mu.txt', 'w');
fp2 = fopen('vmax_mv.txt', 'w');
fp3 = fopen('vmax_mu_satu.txt', 'w');
fp4 = fopen('vmax_mv_satu.txt', 'w');
fp5 = fopen('vmax_m0_satu.txt', 'w');

% arcq_0p1 = 1.29;
% arcq_0p05 = 1.65;
% arcq_0p01 = 2.33;
% arcq_1e-3 = 3.08;

kk=0;
for snr=snr_start:snr_step:snr_end
    kk = kk+1;
	sigma_m(kk) = sqrt((10^(-snr/10.0)/code_rate)/2);
end
N = kk;

vmax_mu = zeros(N, evl_num);
vmax_mv = zeros(N, evl_num);

for ii=1:1:N
    [mu, mv, m0] = evl_ir(sigma_m(ii), lam, rou, evl_num);
    for jj=1:1:evl_num
        vmax_mu(ii, jj) = arcq*sqrt(2*mu(jj))+mu(jj);   
        vmax_mv(ii, jj) = arcq*sqrt(2*mv(jj))+mv(jj);
        vmax_mu_satu(ii, jj) = vmax_mu(ii, jj);
        vmax_mv_satu(ii, jj) = vmax_mv(ii, jj);
        if (vmax_mu(ii, jj) >= satu)
            vmax_mu_satu(ii, jj) = satu;
        end
        if (vmax_mv(ii, jj) >= satu)
            vmax_mv_satu(ii, jj) = satu;
        end
        fprintf(fp1, ' %f', vmax_mu(ii, jj));
        fprintf(fp2, ' %f', vmax_mv(ii, jj));
        fprintf(fp3, ' %f', vmax_mu_satu(ii, jj));
        fprintf(fp4, ' %f', vmax_mv_satu(ii, jj));
    end
    vmax_m0(ii) = arcq*sqrt(2*m0)+m0;
    vmax_m0_satu(ii) = vmax_m0(ii);
    if (vmax_m0(ii) >= satu)
            vmax_m0_satu(ii) = satu;
    end
    fprintf(fp1, '\n');
    fprintf(fp2, '\n');
    fprintf(fp3, '\n');
    fprintf(fp4, '\n');
    fprintf(fp5, '\n');
end
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);

fprintf(fp, 'double vmax_m0_evl[%d]={', N);
for ii=1:1:N-1
    fprintf(fp, '%f, ', vmax_m0(ii));
end
ii=ii+1;
fprintf(fp, '%f};\n', vmax_m0(ii));

fprintf(fp, 'double vmax_m0_satu_evl[%d]={', N);
for ii=1:1:N-1
    fprintf(fp, '%f, ', vmax_m0_satu(ii));
end
ii=ii+1;
fprintf(fp, '%f};\n', vmax_m0_satu(ii));
    
    
fprintf(fp, 'double vmax_mu_evl[%d][%d]={', N, evl_num);
for ii=1:1:N-1
    fprintf(fp, '{');
    for jj=1:1:evl_num-1
        fprintf(fp, '%f, ', vmax_mu(ii, jj));
    end
    jj=jj+1;
    fprintf(fp, '%f},\n', vmax_mu(ii, jj));
end
ii=ii+1;
fprintf(fp, '{');
for jj=1:1:evl_num-1
	fprintf(fp, '%f, ', vmax_mu(ii, jj));
end
jj=jj+1;
fprintf(fp, '%f}\n', vmax_mu(ii, jj));
fprintf(fp, '};\n');

fprintf(fp, 'double vmax_mv_evl[%d][%d]={', N, evl_num);
for ii=1:1:N-1
    fprintf(fp, '{');
    for jj=1:1:evl_num-1
        fprintf(fp, '%f, ', vmax_mv(ii, jj));
    end
    jj=jj+1;
    fprintf(fp, '%f},\n', vmax_mv(ii, jj));
end
ii=ii+1;
fprintf(fp, '{');
for jj=1:1:evl_num-1
	fprintf(fp, '%f, ', vmax_mv(ii, jj));
end
jj=jj+1;
fprintf(fp, '%f}\n', vmax_mv(ii, jj));
fprintf(fp, '};\n');

fprintf(fp, 'double vmax_mu_satu_evl[%d][%d]={', N, evl_num);
for ii=1:1:N-1
    fprintf(fp, '{');
    for jj=1:1:evl_num-1
        fprintf(fp, '%f, ', vmax_mu_satu(ii, jj));
    end
    jj=jj+1;
    fprintf(fp, '%f},\n', vmax_mu_satu(ii, jj));
end
ii=ii+1;
fprintf(fp, '{');
for jj=1:1:evl_num-1
	fprintf(fp, '%f, ', vmax_mu_satu(ii, jj));
end
jj=jj+1;
    fprintf(fp, '%f}\n', vmax_mu_satu(ii, jj));
fprintf(fp, '};\n');

fprintf(fp, 'double vmax_mv_satu_evl[%d][%d]={', N, evl_num);
for ii=1:1:N-1
    fprintf(fp, '{');
    for jj=1:1:evl_num-1
        fprintf(fp, '%f, ', vmax_mv_satu(ii, jj));
    end
    jj=jj+1;
    fprintf(fp, '%f},\n', vmax_mv_satu(ii, jj));
end
ii=ii+1;
fprintf(fp, '{');
for jj=1:1:evl_num-1
	fprintf(fp, '%f, ', vmax_mv_satu(ii, jj));
end
jj=jj+1;
fprintf(fp, '%f}\n', vmax_mv_satu(ii, jj));
fprintf(fp, '};\n');

fclose(fp);

%end
