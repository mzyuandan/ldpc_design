clc;
clear all;
close all;
%%%%%
% H_c_0=[    1   0   0   1   0   0   0   0   0   1   0   1   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
%            1   0   0   0   1   0   1   0   1   0   1   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0;
%            0   1   0   0   0   0   1   0   0   1   0   0   0   1   1   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0;
%            0   0   1   0   1   0   0   0   1   0   1   0   0   1   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0;
%            1   0   0   0   0   1   0   1   0   0   0   0   0   0   1   1   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0;
%            0   1   1   0   1   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0;
%            1   0   0   0   0   1   0   0   1   0   0   1   0   1   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0;
%            0   0   1   0   1   0   0   1   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0;
%            1   0   0   0   0   0   1   0   1   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0;
%            0   0   1   0   1   0   0   0   1   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0;
%            0   1   1   0   0   0   1   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0;
%            0   0   0   1   1   0   0   0   1   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0;
%            1   0   0   0   0   0   0   1   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0;
%            0   0   1   0   1   0   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0;
%            0   0   1   1   0   0   0   0   1   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1;
%            1   0   0   0   0   1   1   0   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1;];
% z_max_0 = 360;


% Dc = [7,7,7,7,7,7,7,7,7,7,7,7,6,7,7,6];

% M = 9;
% N = 27;
% var_vect = [8, 7, 12];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 5, 4];
% ds = [4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2];
% var_vect = [8, 13, 0, 0, 6];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 7, 0, 2];
% ds = [3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2];

% M = 6;
% N = 30;
% var_vect = [5, 12, 13];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2];
% ds = [3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 2, 2, 2, 2, 2];
% var_vect = [5, 17, 0, 0, 8];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 1];
% ds = [3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 2, 2, 2, 2, 2];

% M = 18;
% N = 27;
% % var_vect = [17, 3, 0, 0, 7];
% % chk_vect = [0, 0, 5, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% % ds = [6, 3, 6, 6, 6, 3, 6, 6, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
% % var_vect = [17, 2, 0, 0, 8];
% % chk_vect = [0, 0, 10, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% % ds = [6, 6, 6, 6, 3, 6, 6, 6, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
% % var_vect = [17, 3, 0, 0, 7];
% % chk_vect = [0, 0, 11, 1, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% % ds = [3, 3, 6, 6, 6, 6, 6, 6, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
% var_vect = [17, 6, 0, 0, 4];
% chk_vect = [0, 0, 0, 13, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% ds = [3, 3, 3, 3, 3, 12, 12, 12, 12, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];


% M = 9;
% N = 27;
% var_vect = [8, 7, 12];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 5, 4];
% ds = [4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2];
% var_vect = [8, 13, 0, 0, 6];
% chk_vect = [0, 0, 0, 0, 0, 0, 0, 7, 0, 2];
% ds = [3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2];

M = 13;
N = 18;
var_vect = [12, 1, 0, 0, 5];
chk_vect = [0, 4, 0, 9];
ds = [6, 6, 6, 6, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];


K = N - M;
H_b = zeros(M,N);
ds_m = length(var_vect) + 1;
dc_m = length(chk_vect) + 1;
dc_ini = 2 * ones(1, M);
dc_ini(ceil(M/2)) = 3;


yps_tp = find(chk_vect ~= 0) + 1;
yps_set = [yps_tp(1), yps_tp(1)+1, yps_tp(1)+2];
xps_set = zeros(1, length(yps_set));
for ii=1:1:length(yps_set)
    xps_set(ii) = chk_vect(yps_set(ii)-1);
end
chk_vect_g = zeros(1, length(chk_vect));
chk_vect_g(1) = M - 1;
chk_vect_g(2) = 1;
chk_max = chk_vect(yps_set(2)-1) + chk_vect(yps_set(3)-1);


%%校验部分赋初值
for ii=1:1:M-1
    for jj=ii+K:1:ii+K+1
        H_b(ii,jj) = 1;
    end
end
H_b(M,N) = 1;
H_b(ceil(M/2),K+1) = 1;
H_b(M,K+1) = 1;
%%变量节点度分布
% ds = [7,3,7,3,7,3,7,3,7,3,7,3,3,7,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];


flg = 0;

while (flg==0)
    %%开始PEG算法
    
    while(1)
        %%校验节点初始度分布
        dc = dc_ini;
        yps_min_rc = 0;
        for ii=1:1:M
            for jj=1:1:K
                H_b(ii,jj)=0;
            end
        end
        rand_jj = zeros(1,K);
        for jj=1:1:K
%             rjj_flg = 0;
%             while (rjj_flg==0)
%                 jjx = randi(K, 1);
%                 if (rand_jj(jjx)==0)
%                     rand_jj(jjx)=1;
%                     rjj_flg = 1;
%                 end
%             end
            jjx =jj;
            for kk=1:1:ds(jjx)
                if (kk == 1)
                    dc_min = min(dc);
                    
                    if ((dc_min==yps_set(1)) && ((chk_vect_g(yps_set(2)-1)+chk_vect_g(yps_set(3)-1))==chk_max))
                        dc_min = dc_min + 1;
                        yps_min_rc = 1;
                    end
                    
                    
                    k1 = find(dc==dc_min);       %%第一条边，寻找度最小的校验节点
                    
                    chk_flg = 0;
                    while (chk_flg==0)
                        xtp = randi(length(k1), 1);
                        tp = dc(k1(xtp));
                        chk_flg = 1;
                        if (dc_min<=yps_set(1))
                            if ((tp+1>yps_set(1)) && ((chk_vect_g(yps_set(2)-1)+chk_vect_g(yps_set(3)-1))==chk_max))
                                chk_flg = 0;
                            end
                        else
                            if ((tp==yps_set(2)) && (chk_vect_g(yps_set(3)-1)==chk_vect(yps_set(3)-1)))
                                chk_flg = 0;
                            end
                        end
                    end
                    chk_vect_g(tp) = chk_vect_g(tp) + 1;
                    chk_vect_g(tp-1) = chk_vect_g(tp-1) - 1;
                    
                    H_b(k1(xtp),jjx) = 1;
                    dc(k1(xtp)) = dc(k1(xtp))+1;
                else
                    ll = 1;                         %%深度为1   
                    flag = zeros(M,M);              %%初存每一层与该变量节点关联的校验节点，flag=1表示有关联，flag=0表示无关联
                    row = find(H_b(:,jjx));
                    flag(ll,row) = 1;
                    while(1)
                        ll=ll+1;
                        H_row = H_b(row,:);
                        if (length(row)==1)
                            col = find(H_row);                    %%行搜索
                        else
                            col = find(any(H_row(:,1:N)));
                        end

                        H_col = (H_b(:,col))';
                        if (length(col)==1)
                            row = find(H_col);                     %%列搜索
                        else
                            row = find(any(H_col(:,1:M)));
                        end
                        flag(ll,row) = 1;                           %%与该变量节点关联的校验节点

                        if (ll>=2)
                            %%第一种情况，关联的校验节点不再增加
                            if (sum(flag(ll-1,:))==sum(flag(ll,:)) && sum(flag(ll-1,:))<M)
                                flag_0_all = find(flag(ll,:)==0);
                                dc_0 = dc(flag_0_all);
                                flag_0 = find(dc_0==min(dc_0));      %%寻找度最小的校验节点
                                m = length(flag_0);
        %                         while(1)
                                    
                                    
                                    chk_flg = 0;
                                    while (chk_flg==0)
                                        n = randi(m, 1, 1); %randint(1,1,m)+1
                                        tp = dc(flag_0_all(flag_0(n)));
                                        chk_flg = 1;
                                        if ((tp==yps_set(1)) && ((chk_vect_g(yps_set(2)-1)+chk_vect_g(yps_set(3)-1))==chk_max))
                                            chk_flg = 0;
                                        end
                                    end
                                    chk_vect_g(tp) = chk_vect_g(tp) + 1;
                                    chk_vect_g(tp-1) = chk_vect_g(tp-1) - 1;
                                    
                                    H_b(flag_0_all(flag_0(n)),jjx) = 1;
                                    dc(flag_0_all(flag_0(n))) = dc(flag_0_all(flag_0(n)))+1;
        %                             if (dc(flag_0_all(flag_0(n))) <= 7)
        %                                 break;
        %                             end
        %                         end
                                break;
                            %%第二种情况，校验节点饱和
                            elseif (sum(flag(ll-1,:))<M && sum(flag(ll,:))==M)
                                flag_0_all = find((flag(ll-1,:) & flag(ll,:))==0);
                                dc_0 = dc(flag_0_all);
                                flag_0 = find(dc_0==min(dc_0));                 %%寻找度最小的校验节点
                                m = length(flag_0);
                                    
                                chk_flg = 0;
                                tested = zeros(1, m);
                                dcpp = 0;
                                while (chk_flg==0)
                                    n = randi(m, 1, 1); %randint(1,1,m)+1
                                    tp = dc(flag_0_all(flag_0(n)));
                                    chk_flg = 1;
                                    if (dcpp==0)
                                        if ((tp==yps_set(1)) && ((chk_vect_g(yps_set(2)-1)+chk_vect_g(yps_set(3)-1))==chk_max))
                                            chk_flg = 0;
                                        end
                                    else
                                        if (dcpp==1)
                                            if ((tp==yps_set(2)) && (chk_vect_g(yps_set(3)-1)==chk_vect(yps_set(3)-1)))
                                                chk_flg = 0;
                                            end
                                        end
                                    end
                                   
                                    tested(n) = 1;
                                    if (sum(tested)==m && chk_flg==0)
                                        dcpp = dcpp + 1;
                                        flag_0 = find(dc_0==min(dc_0)+dcpp);                 %%寻找度最小的校验节点
                                        m = length(flag_0);

                                        chk_flg = 0;
                                        tested = zeros(1, m);
                                        
                                    end

                                end
                                chk_vect_g(tp) = chk_vect_g(tp) + 1;
                                chk_vect_g(tp-1) = chk_vect_g(tp-1) - 1;

                                H_b(flag_0_all(flag_0(n)),jjx) = 1;
                                dc(flag_0_all(flag_0(n))) = dc(flag_0_all(flag_0(n)))+1;
        %                              if (dc(flag_0_all(flag_0(n))) <= 7)
        %                                  break;
        %                              end
        %                         end
                                break;
                            end
                        end
                    end
                end
            end
        end

        dpm = dc(1);
        dpl = 0;
        flag_dc = 1;
        for nn=2:1:M
            dpf = 0;
            for mm=1:1:dpl
                if (dc(2)==dpm(dpl))
                    dpf = 1;
                    break;
                end
            end
            if (dpf==0)
                dpl = dpl + 1;
                dpm = [dpm, dc(nn)];
            end
            if (dpl>2)
                flag_dc = 0; 
                break;
            end
        end

        if (flag_dc == 1)
            break;
        end
    end
                        
    for ii=1:1:M
        for jj=1:1:N
            fprintf('%4d',H_b(ii,jj));
        end
        fprintf('\n');
    end
    fprintf('\n');
    % for ii=1:1:M
    fprintf('%4d',dc);
    fprintf('\n');   

    chk_vec_t = zeros(1, length(chk_vect));
    yps1 = min(dc);
    yps2 = yps1 + 1;
    yps3 = yps1 + 2;
    xps1 = sum(dc==yps1);
    xps2 = sum(dcds==yps2);
    xps3 = sum(dc==yps3);
    chk_vec_t(yps1-1) = xps1;
    chk_vec_t(yps2-1) = xps2;
    chk_vec_t(yps3-1) = xps3;

    if (sum(chk_vec_t==chk_vect)==length(chk_vect))
        flg = 1;
    end
end

aa = 0;
                
                
            
            
           
        