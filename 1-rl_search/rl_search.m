clc;
clear all;
close all;


mb = 13;
nb = 18;
th_st1 = 1.28;
th_st2 = 1.28;


kb = nb - mb;
dv = mb;
dc = nb;
evl_num = 50;
Qval_t = 10;
CodeRate = kb / nb;


th_max = 0;
cnt = 0;
% for ii=1:1:kb+1
%     dv = 4;
%     var_lin = 2:1:dv;
%     var_lin = var_lin';
%     x1 = ii;
%     var_vect = zeros(1, dv-1);
%     var_vect(1) = mb-1;
%     var_vect(2) = x1;
%     var_vect(dv-1) = kb+1-x1;
% %     if (sum(var_vect)~=nb)
% %         break;
% %     end
%     
%     weight = var_vect * var_lin;
%     
%     for jj=0:1:mb
%         x21 = jj;
%         xx2 = mb - x21;
%         for kk=0:1:xx2
%             x22 = kk;
%             y2 = (weight + x21 - x22 - mb) / mb;
%             y2_int = floor(y2);
%             if (y2==y2_int)    
%                 dc = y2+2;
%                 chk_lin = 2:1:dc;
%                 chk_lin = chk_lin';
%                 chk_vect = zeros(1, dc-1);
%                 chk_vect(y2-1) = x21;
%                 chk_vect(y2) = mb - x21 - x22;
%                 chk_vect(y2+1) = x22;
% %                 if (sum(chk_vect)~=mb)
% %                     break;
% %                 end
% %                 rou = roux / mb;
% %                 weight1 = chk_vect * chk_lin;
% %                 if (weight1~=weight)
% %                     break;
% %                 end
%                 
%                 lam = zeros(1, dv-1);
%                 for ll=2:1:dv
%                     lam(ll-1) = ll*var_vect(ll-1) / weight;
%                 end
% 
%                 rou = zeros(1, dc-1);
%                 for ll=2:1:dc
%                     rou(ll-1) = ll*chk_vect(ll-1) / weight;
%                 end
% 
% 
%                 threshold = threshold_cal_ir_new1( lam, rou, evl_num, Qval_t, 0.0001, th_st1 );
%                 fprintf('dv=%d, x1=%d, x21=%d, x22=%d, y2=%d, th=%.4f\n', dv, x1, x21, x22, y2, threshold);
%                 if (threshold>th_max)
%                     th_max = threshold;
%                     lamt = lam;
%                     rout = rou;
%                     var_vect_t = var_vect;
%                     chk_vect_t = chk_vect;
%                     fprintf('th_max = %.4f;\n', th_max);
%                 end
%             end
%         end
%     end
%     cnt = cnt + 1;
% end

if (mb>=6)
    for ii=1:1:kb+1
        dv = 6;
        var_lin = 2:1:dv;
        var_lin = var_lin';
        x1 = ii;
        var_vect = zeros(1, dv-1);
        var_vect(1) = mb-1;
        var_vect(2) = x1;
        var_vect(dv-1) = kb+1-x1;
%         if (sum(var_vect)~=nb)
%             break;
%         end

        weight = var_vect * var_lin;

        for jj=0:1:mb
            x21 = jj;
            xx2 = mb - x21;
            for kk=0:1:xx2
                x22 = kk;
                y2 = (weight + x21 - x22 - mb) / mb;
                y2_int = floor(y2);
                if (y2==y2_int)    
                    dc = y2+2;
                    chk_lin = 2:1:dc;
                    chk_lin = chk_lin';
                    chk_vect = zeros(1, dc-1);
                    chk_vect(y2-1) = x21;
                    chk_vect(y2) = mb - x21 - x22;
                    chk_vect(y2+1) = x22;
%                     if (sum(chk_vect)~=mb)
%                         break;
%                     end
    %                 rou = roux / mb;
%                     weight1 = chk_vect * chk_lin;
%                     if (weight1~=weight)
%                         break;
%                     end

                    lam = zeros(1, dv-1);
                    for ll=2:1:dv
                        lam(ll-1) = ll*var_vect(ll-1) / weight;
                    end

                    rou = zeros(1, dc-1);
                    for ll=2:1:dc
                        rou(ll-1) = ll*chk_vect(ll-1) / weight;
                    end


                    threshold = threshold_cal_ir_new1( lam, rou, evl_num, Qval_t, 0.0001, th_st2 );
                    fprintf('dv=%d, x1=%d, x21=%d, x22=%d, y2=%d, th=%.4f\n', dv, x1, x21, x22, y2, threshold);
                    if (threshold>th_max)
                        th_max = threshold;
                        lamt = lam;
                        rout = rou;
                        var_vect_t = var_vect;
                        chk_vect_t = chk_vect;
                        fprintf('th_max = %.4f;\n', th_max);
                    end
                end
            end
        end
        cnt = cnt + 1;
    end
end

DateString = datestr(datetime('now'));
fid = fopen('result.txt', 'a');


fprintf(fid, '----------------%s----------------\n', DateString);
fprintf(fid, 'mb=%d nb=%d th=%.4f\n', mb, nb, th_max);
fprintf(fid, 'var_wgt:\t');
fprintf(fid, '%d\t', var_lin);
fprintf(fid, '\n');
fprintf(fid, 'var_cnt:\t');
fprintf(fid, '%d\t', var_vect_t);
fprintf(fid, '\n');
fprintf(fid, 'chk_wgt:\t');
fprintf(fid, '%d\t', chk_lin);
fprintf(fid, '\n');
fprintf(fid, 'chk_cnt:\t');
fprintf(fid, '%d\t', chk_vect_t);
fprintf(fid, '\n\n\n');

aa=0;



