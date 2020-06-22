%%%%基于16*32基矩阵生成Z=360,180,90,45都不含4,6环的矩阵
clc;
clear all;
close all;
%%%%%新生成的码字基矩阵  0对应矩阵中的-1,1对应非零值
QC_H_ori = dlmread('hb_9x27A1.txt');
[M, N] = size(QC_H_ori);
z_max_s = 512;


K = N - M;
CodeRate = K / N;
xpos_s = ceil(M / 2);

%矩阵中0替换为-1，不需要则屏蔽掉此段程序
for ii=1:1:M
    for jj=1:1:N
        if (QC_H_ori(ii,jj)==0)
            QC_H_ori(ii,jj)=-1;
        else
            QC_H_ori(ii,jj) = QC_H_ori(ii,jj);
        end
            fprintf('%4d',QC_H_ori(ii,jj));
     end
    fprintf('\n');
end
fprintf('\n');

%%
H_c = QC_H_ori;
z_max = z_max_s;
xpos = xpos_s;
Z1 = 60;
Z2 = 40;
Z3 = 20;
Z4 = 10;
g4 = 0;
g6 = 0;
g4_pos = zeros(1000, 8);
g6_pos = zeros(3000,12);
[M, N] = size(H_c);
K = N - M;
% H_c = zeros(M,N);                      %%%输出的Zmax=256的矩阵
H_c_b0 = -1.*ones(M,N); 
H_c_b1 = -1.*ones(M,N);                  %%Z1=360的准循环移位矩阵
H_c_b2 = -1.*ones(M,N);                  %%Z2=180的准循环移位矩阵
H_c_b3 = -1.*ones(M,N);                  %%Z3=90的准循环移位矩阵
H_c_b4 = -1.*ones(M,N);                  %%Z4=45的准循环移位矩阵
zz = zeros(M,1);                         %%存储校验节点的度

g4s = 1000*ones(1,4);
g6s = 1000*ones(1,4);
set_g4 = [0 0 0 0];
set_g6 = [0 0 0 80];


  %%矩阵中0变成-1
for ii=1:1:M
    for jj=1:1:N
        if (H_c(ii,jj)==0)
            H_c(ii,jj)=-1;
        else
            H_c(ii,jj) = H_c(ii,jj);
        end
            fprintf('%4d',H_c(ii,jj));
     end
    fprintf('\n');
end
fprintf('\n');

  %%检测基矩阵的4.6环
[g4,g6,g4_pos,g6_pos,xx1, xx2] = my_girth_4_6(H_c,z_max,z_max,g4,g6,g4_pos,g6_pos, 1);


%% 赋初值，输出无四环
  H_c_b0(1,K+2) = 0;
  H_c_b1(1,K+2) = 0;
  H_c_b2(1,K+2) = 0;
  H_c_b3(1,K+2) = 0;
  H_c_b4(1,K+2) = 0;
  for ii=2:1:M-1
      for jj = ii+K:1:ii+K+1
          H_c_b0(ii,jj) = 0;
          H_c_b1(ii,jj) = 0;
          H_c_b2(ii,jj) = 0;
          H_c_b3(ii,jj) = 0;
          H_c_b4(ii,jj) = 0;
      end
  end
  H_c_b0(M,N) = 0;
  H_c_b1(M,N) = 0;
  H_c_b2(M,N) = 0;
  H_c_b3(M,N) = 0;
  H_c_b4(M,N) = 0;
  
  
  
 while(1)
     if (girth_cmp(g4s, set_g4)==1)
         break;                         %%无4环，输出
     end
     if (girth_cmp(g4s, set_g4)==0)
        gg = fix(z_max*rand);
        H_c_b0(1,K+1) = gg;
        H_c_b1(1,K+1) = floor(gg*Z1/z_max);
        H_c_b2(1,K+1) = floor(gg*Z2/z_max);
        H_c_b3(1,K+1) = floor(gg*Z3/z_max);
        H_c_b4(1,K+1) = floor(gg*Z4/z_max);
        H_c_b0(M,K+1) = gg;
        H_c_b1(M,K+1) = floor(gg*Z1/z_max);
        H_c_b2(M,K+1) = floor(gg*Z2/z_max);
        H_c_b3(M,K+1) = floor(gg*Z3/z_max);
        H_c_b4(M,K+1) = floor(gg*Z4/z_max);
        
        H_c_b0(xpos,K+1) = fix(z_max*rand);
        H_c_b1(xpos,K+1) = floor(H_c_b0(xpos,17)*Z1/z_max);
        H_c_b2(xpos,K+1) = floor(H_c_b0(xpos,17)*Z2/z_max);
        H_c_b3(xpos,K+1) = floor(H_c_b0(xpos,17)*Z3/z_max);
        H_c_b4(xpos,K+1) = floor(H_c_b0(xpos,17)*Z4/z_max);
        
        for ii=1:1:g4
            cnt = 0;
            while(1)
                cnt = cnt + 1;
                if (cnt>=1000)
                    break;
                end
                for jj=1:2:7
                    if (g4_pos(ii,jj+1)>=K+1)
                        continue;
                    else 
                        H_c_b0(g4_pos(ii,jj),g4_pos(ii,jj+1)) = fix(z_max*rand);
                        H_c_b1(g4_pos(ii,jj),g4_pos(ii,jj+1)) = floor(H_c_b0(g4_pos(ii,jj),g4_pos(ii,jj+1))*Z1/z_max);
                        H_c_b2(g4_pos(ii,jj),g4_pos(ii,jj+1)) = floor(H_c_b0(g4_pos(ii,jj),g4_pos(ii,jj+1))*Z2/z_max);
                        H_c_b3(g4_pos(ii,jj),g4_pos(ii,jj+1)) = floor(H_c_b0(g4_pos(ii,jj),g4_pos(ii,jj+1))*Z3/z_max);
                        H_c_b4(g4_pos(ii,jj),g4_pos(ii,jj+1)) = floor(H_c_b0(g4_pos(ii,jj),g4_pos(ii,jj+1))*Z4/z_max);
                    end
                end                      
                if ((mod((H_c_b1(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b1(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b1(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b1(g4_pos(ii,7),g4_pos(ii,8))),Z1)==0) || ...
                    (mod((H_c_b2(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b2(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b2(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b2(g4_pos(ii,7),g4_pos(ii,8))),Z2)==0) || ...
                    (mod((H_c_b3(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b3(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b3(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b3(g4_pos(ii,7),g4_pos(ii,8))),Z3)==0) || ...
                    (mod((H_c_b4(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b4(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b4(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b4(g4_pos(ii,7),g4_pos(ii,8))),Z4)==0));
                    continue;
                else
                    break;
                end
            end
        end
        
  %%对不在环中的未被替换的"1"元素赋值     
        for ii=1:1:M
            for jj=1:1:K
                for pp=1:1:g4
                    for qq=1:2:7
                        flag_4_2=0;
                        if (ii==g4_pos(pp,qq) && jj==g4_pos(pp,qq+1))      %%除去已赋过的1元素
                            flag_4_2=1;
                            break;
                        end
                    end
                    if (flag_4_2==1)
                        break;
                    end
                end
                if(flag_4_2==1)
                   continue;
                else
                    if (H_c(ii,jj)>=0)
                       H_c_b0(ii,jj) = fix(z_max*rand);
                       H_c_b1(ii,jj) = floor(H_c_b0(ii,jj)*Z1/z_max);
                       H_c_b2(ii,jj) = floor(H_c_b0(ii,jj)*Z2/z_max);
                       H_c_b3(ii,jj) = floor(H_c_b0(ii,jj)*Z3/z_max);
                       H_c_b4(ii,jj) = floor(H_c_b0(ii,jj)*Z4/z_max);
                     end
                end
            end
        end
     end
     [g4,g6,g4_pos,g6_pos,g4s(1),g6s(1)] = my_girth_4_6(H_c_b1,Z1,z_max,g4,g6,g4_pos,g6_pos, 1);
     [g4,g6,g4_pos,g6_pos,g4s(2),g6s(2)] = my_girth_4_6(H_c_b2,Z2,z_max,g4,g6,g4_pos,g6_pos, 0);
     [g4,g6,g4_pos,g6_pos,g4s(3),g6s(3)] = my_girth_4_6(H_c_b3,Z3,z_max,g4,g6,g4_pos,g6_pos, 0);
     [g4,g6,g4_pos,g6_pos,g4s(4),g6s(4)] = my_girth_4_6(H_c_b4,Z4,z_max,g4,g6,g4_pos,g6_pos, 0);
 end

for ii=1:1:M
    for jj=1:1:N
        fprintf('%4d,', H_c_b0(ii,jj));
    end
    fprintf('\n');
end
% fprintf('\n');
% for ii=1:1:M
%     for jj=1:1:N
%         fprintf('%4d,', H_c_b1(ii,jj));
%     end
%     fprintf('\n');
% end
% fprintf('\n');
% for ii=1:1:M
%     for jj=1:1:N
%         fprintf('%4d,', H_c_b2(ii,jj));
%     end
%     fprintf('\n');
% end
% fprintf('\n');
% for ii=1:1:M
%     for jj=1:1:N
%         fprintf('%4d,', H_c_b3(ii,jj));
%     end
%     fprintf('\n');
% end
% fprintf('\n');
% for ii=1:1:M
%     for jj=1:1:N
%         fprintf('%4d,', H_c_b4(ii,jj));
%     end
%     fprintf('\n');
% end


%% 消4,6环
g4_min = 1000000;
g6_min = 1000000;
cnt_g6 = 0;
cnt_g4 = 0;
g6_flg = 0;

 while (1)
    if (girth_cmp(g4s, set_g4)==1 && girth_cmp(g6s, set_g6)==1)                        %%%无4,6环，输出结果
        for ii=1:1:M
            for jj=1:1:N
               fprintf('%4d',H_c_b0(ii,jj));
            end
            fprintf('\n');
        end
        fprintf('\n');
        break;
    end
    
    
    
    %%消6环
    while(1)
        if (girth_cmp(g4s, set_g4)==0)
            break;                         %%有4环，先消4环
        end
        if (girth_cmp(g6s, set_g6)==1)                         %%无6环，跳出
            break;
        end
        if (girth_cmp(g6s, set_g6)==0)  
            g6_flg = 1;
            for ii=1:1:M
                for jj=1:1:N
                    if (H_c_b1(ii,jj)>=0)
                        zz(ii)=zz(ii)+1;              %%计算校验节点的度
                    end
                end
            end
        
                            
            for ii=1:1:g6
               gogo = (mod((H_c_b1(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b1(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b1(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b1(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b1(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b1(g6_pos(ii,11),g6_pos(ii,12))),Z1)==0) || ...
                          (mod((H_c_b2(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b2(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b2(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b2(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b2(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b2(g6_pos(ii,11),g6_pos(ii,12))),Z2)==0) || ...
                          (mod((H_c_b3(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b3(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b3(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b3(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b3(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b3(g6_pos(ii,11),g6_pos(ii,12))),Z3)==0) || ...
                          (mod((H_c_b4(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b4(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b4(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b4(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b4(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b4(g6_pos(ii,11),g6_pos(ii,12))),Z4)==0);
               while(gogo)
                      if (zz(g6_pos(ii,1))>zz(g6_pos(ii,5)))
                          if (zz(g6_pos(ii,5))>zz(g6_pos(ii,9)))                %%改变校验节点度小的Z
                              ll = 9;
                              oo = 5;
                          else
                              ll = 5;
                              oo = 9;
                          end
                      else
                          if (zz(g6_pos(ii,1))>zz(g6_pos(ii,9)))
                              ll = 9;
                              oo = 1;
                          else
                              ll = 1;
                              oo = 9;
                          end
                      end
                      if (g6_pos(ii,ll+1)>K)
                           ll = ll + 2;
                      else if (g6_pos(ii,ll+1)>K)
                               ll = oo;
                          end
                      end
                      

                      H_c_b0(g6_pos(ii,ll),g6_pos(ii,ll+1)) = fix(z_max*rand);
                      H_c_b1(g6_pos(ii,ll),g6_pos(ii,ll+1)) = floor(H_c_b0(g6_pos(ii,ll),g6_pos(ii,ll+1))*Z1/z_max);
                      H_c_b2(g6_pos(ii,ll),g6_pos(ii,ll+1)) = floor(H_c_b0(g6_pos(ii,ll),g6_pos(ii,ll+1))*Z2/z_max);
                      H_c_b3(g6_pos(ii,ll),g6_pos(ii,ll+1)) = floor(H_c_b0(g6_pos(ii,ll),g6_pos(ii,ll+1))*Z3/z_max);
                      H_c_b4(g6_pos(ii,ll),g6_pos(ii,ll+1)) = floor(H_c_b0(g6_pos(ii,ll),g6_pos(ii,ll+1))*Z4/z_max);
%                       if (H_c_b(1,17)~=gg)
%                            
%                       end
        
                      if ((mod((H_c_b1(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b1(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b1(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b1(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b1(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b1(g6_pos(ii,11),g6_pos(ii,12))),Z1)==0) || ...
                          (mod((H_c_b2(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b2(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b2(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b2(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b2(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b2(g6_pos(ii,11),g6_pos(ii,12))),Z2)==0) || ...
                          (mod((H_c_b3(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b3(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b3(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b3(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b3(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b3(g6_pos(ii,11),g6_pos(ii,12))),Z3)==0) || ...
                          (mod((H_c_b4(g6_pos(ii,1),g6_pos(ii,2)) - H_c_b4(g6_pos(ii,3),g6_pos(ii,4)) + H_c_b4(g6_pos(ii,5),g6_pos(ii,6)) -  ...
                               H_c_b4(g6_pos(ii,7),g6_pos(ii,8)) + H_c_b4(g6_pos(ii,9),g6_pos(ii,10)) - H_c_b4(g6_pos(ii,11),g6_pos(ii,12))),Z4)==0));
                          continue;
                      else
                          break;
                      end
               end
            end
         end
         [g4,g6,g4_pos,g6_pos,g4s(1),g6s(1)] = my_girth_4_6(H_c_b1,Z1,z_max,g4,g6,g4_pos,g6_pos, 1);
         [g4,g6,g4_pos,g6_pos,g4s(2),g6s(2)] = my_girth_4_6(H_c_b2,Z2,z_max,g4,g6,g4_pos,g6_pos, 0);
         [g4,g6,g4_pos,g6_pos,g4s(3),g6s(3)] = my_girth_4_6(H_c_b3,Z3,z_max,g4,g6,g4_pos,g6_pos, 0);
         [g4,g6,g4_pos,g6_pos,g4s(4),g6s(4)] = my_girth_4_6(H_c_b4,Z4,z_max,g4,g6,g4_pos,g6_pos, 0);
     
         if(g4s(4)==0 && g6s(4)<g6_min)
             g6_min = g6s(4);
             H_c_rec = H_c_b0;
         end
         cnt_g6 = cnt_g6 + 1;
    end
   
    %%消4环
    while(1)
        if (girth_cmp(g4s, set_g4)==1)                   %%无4环，跳出
            break;
        end
        if (girth_cmp(g4s, set_g4)==0)                      
            for ii=1:1:M
                for jj=1:1:N
                    if (H_c_b1(ii,jj)>=0)
                        zz(ii)=zz(ii)+1;     %%计算校验节点的度     
                    end
                end
            end                     
            for ii=1:1:g4
               while(1)
                      if (zz(g4_pos(ii,1)) > zz(g4_pos(ii,5)))               %%%%改变校验节点度小的Z
                          ll = 5;
                      else
                          ll = 1;
                      end
                      if(g4_pos(ii,ll+1)>K)
                          ll=ll+2;
                      end
                          
                      H_c_b0(g4_pos(ii,ll),g4_pos(ii,ll+1)) = fix(z_max*rand);
                      H_c_b1(g4_pos(ii,ll),g4_pos(ii,ll+1)) = floor(H_c_b0(g4_pos(ii,ll),g4_pos(ii,ll+1))*Z1/z_max);
                      H_c_b2(g4_pos(ii,ll),g4_pos(ii,ll+1)) = floor(H_c_b0(g4_pos(ii,ll),g4_pos(ii,ll+1))*Z2/z_max);
                      H_c_b3(g4_pos(ii,ll),g4_pos(ii,ll+1)) = floor(H_c_b0(g4_pos(ii,ll),g4_pos(ii,ll+1))*Z3/z_max);
                      H_c_b4(g4_pos(ii,ll),g4_pos(ii,ll+1)) = floor(H_c_b0(g4_pos(ii,ll),g4_pos(ii,ll+1))*Z4/z_max);
                      
                      if ((mod((H_c_b1(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b1(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b1(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b1(g4_pos(ii,7),g4_pos(ii,8))),Z1)==0) || ...
                          (mod((H_c_b2(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b2(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b2(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b2(g4_pos(ii,7),g4_pos(ii,8))),Z2)==0) || ...
                          (mod((H_c_b3(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b3(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b3(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b3(g4_pos(ii,7),g4_pos(ii,8))),Z3)==0) || ...
                          (mod((H_c_b4(g4_pos(ii,1),g4_pos(ii,2)) - H_c_b4(g4_pos(ii,3),g4_pos(ii,4)) + H_c_b4(g4_pos(ii,5),g4_pos(ii,6)) - H_c_b4(g4_pos(ii,7),g4_pos(ii,8))),Z4)==0));
                          continue;
                      else
                          break;
                      end
               end
            end
        end
       [g4,g6,g4_pos,g6_pos,g4s(1),g6s(1)] = my_girth_4_6(H_c_b1,Z1,z_max,g4,g6,g4_pos,g6_pos, 1);
       [g4,g6,g4_pos,g6_pos,g4s(2),g6s(2)] = my_girth_4_6(H_c_b2,Z2,z_max,g4,g6,g4_pos,g6_pos, 0);
       [g4,g6,g4_pos,g6_pos,g4s(3),g6s(3)] = my_girth_4_6(H_c_b3,Z3,z_max,g4,g6,g4_pos,g6_pos, 0);
       [g4,g6,g4_pos,g6_pos,g4s(4),g6s(4)] = my_girth_4_6(H_c_b4,Z4,z_max,g4,g6,g4_pos,g6_pos, 0);
       if (g6_flg==0)
           if (g4s(4)<g4_min)
               g4_min = g4s(4);
               H_c_rec = H_c_b0;
           end
           cnt_g4 = cnt_g4 + 1;
       end
      
    end
 end


    
        
        

            
            
 
          
        
       
    
    





