clc;
clear all;
close all;

mb = 16;
nb = 24;
kb = nb - mb;
evl_num = 50;
BER = 10^(-5);
Qval_t = 2*(qfuncinv(BER)^2);
test_num = 1;
CodeRate = 0.5;

QC_H_ori = dlmread('3GPP_1d3.txt');

threshold = zeros(1, test_num);
SNR_t = zeros(1, test_num);
pos = zeros(test_num, kb);
i1 = 3;
i2 = 8;
i3 = 9;
i4 = 10;
i5 = 11;
i6 = 12;

% QC_H = [QC_H_ori(:,i1) QC_H_ori(:,i2) QC_H_ori(:,i3)...
%     QC_H_ori(:,i4) QC_H_ori(:,i5) QC_H_ori(:,i6) QC_H_ori(:,13:24)];
QC_H = QC_H_ori;
                        
pos = [i1 i2 i3 i4 i5 i6];
col_w = zeros(1, nb);
for jj=1:1:nb
	for ii=1:1:mb
        if (QC_H(ii,jj)~=-1)
            col_w(jj) = col_w(jj) +1;
        end
    end
end
dv = max(col_w);
						
row_w = zeros(1, mb);
for ii=1:1:mb
	for jj=1:1:nb
        if (QC_H(ii,jj)~=-1)
            row_w(ii) = row_w(ii) +1;
        end
	end
end
dc = max(row_w);
						
var_vect = zeros(1, dv-1);
for jj=2:1:dv
    for xx=1:1:nb
        if (col_w(xx)==jj)
            var_vect(jj-1) = var_vect(jj-1) + 1; 
        end
    end
end

chk_vect = zeros(1, dc-1);
for ii=2:1:dc
    for xx=1:1:mb
        if (row_w(xx)==ii)
            chk_vect(ii-1) = chk_vect(ii-1) + 1; 
        end
    end
end

edg_num = 0;
for jj=2:1:dv
    edg_num = edg_num + jj*var_vect(jj-1);
end

lam = zeros(1, dv-1);
for jj=2:1:dv
    lam(jj-1) = jj*var_vect(jj-1) / edg_num;
            end

            rou = zeros(1, dc-1);
for ii=2:1:dc
    rou(ii-1) = ii*chk_vect(ii-1) / edg_num;
end

snr = 0.0:0.2:3.0;
Es_N0 = 10.^snr./10 .* CodeRate;
N0 = 1./(Es_N0);
sigma = N0 ./ 2;

len = length(snr);
mv_all = zeros(1,len);
for ii=1:1:len
    [mv_all(ii), mu, mv, m0] = evl_ir_new(sigma(ii), lam, rou, evl_num);
end
    
    
    

threshold = threshold_cal_ir_new( lam, rou, evl_num, Qval_t );	
SNR_t = 10*log10(1/(2*(threshold^2)));

a=0;



