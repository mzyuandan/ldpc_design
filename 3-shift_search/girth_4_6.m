function [g4,g6,g4_pos,g6_pos,g4s,g6s] = girth_4_6(H_c,Z,g4,g6,g4_pos,g6_pos, rst)
% clc;
% clear all;
% close all;

[M, N] = size(H_c);
H_c_0=zeros(M,N);
for ii=1:1:M
    for jj=1:1:N
        if (H_c(ii,jj)>=0)
            H_c_0(ii,jj)=1;
        else
            H_c_0(ii,jj)=H_c(ii,jj);
        end
%         fprintf('%4d',H_c_0(ii,jj));
    end
%     fprintf('\n');
end

g4s = 0;
g6s = 0;

if (rst)   %%%每次重新计算时清零，但是4种码长4,6环个数要一起累加
    flag_g4 = 0;
    flag_g6 = 0;
    g4 = 0;
    g6 = 0;
    g4_pos = zeros(1000, 8);
    g6_pos = zeros(3000,12);
end

for ii=1:1:M-1
%     fprintf('%d\n', ii);
    one_num1=0;
    one_pos1=zeros(1,10);
    for pp=1:1:N
        if(H_c_0(ii,pp)==1)
            one_num1=one_num1+1;
            one_pos1(one_num1)=pp;
        end
    end
    for jj=ii+1:1:M
        one_num2=0;
        one_pos2=zeros(1,10);
        for qq=1:1:N
            if(H_c_0(jj,qq)==1)
                one_num2=one_num2+1;
                one_pos2(one_num2)=qq;
            end
        end
        
        ovlp12 = 0;
        ovlp12_pos = zeros(1, 8);
        for pp=1:1:one_num1
             for qq=1:1:one_num2
                 if (one_pos1(pp)==one_pos2(qq))
                     ovlp12 = ovlp12 + 1;
                     ovlp12_pos(ovlp12) = one_pos1(pp);
                 end
             end
        end
        if (ovlp12>=2)
%            flag_g4 = 1;
            for pp = 1:1:ovlp12-1
                for qq = pp+1:1:ovlp12
                    if (mod((H_c(ii,ovlp12_pos(pp))-H_c(ii,ovlp12_pos(qq))+H_c(jj,ovlp12_pos(qq))-H_c(jj,ovlp12_pos(pp))),Z)==0)
                       g4 = g4 + 1;
                       g4s = g4s + 1;
                       g4_pos(g4, :) = [ii, ovlp12_pos(pp), ii, ovlp12_pos(qq), jj, ovlp12_pos(qq), jj, ovlp12_pos(pp)];
%                        fprintf('%d %d\n', ii, jj);
%                        for tt=1:1:one_num1
%                           fprintf('%d ', one_pos1(tt));
%                        end
%                        fprintf('\n');
%                        for tt=1:1:one_num2
%                           fprintf('%d ', one_pos2(tt));
%                        end
%                        fprintf('\n');
                    end
                end
            end
        end
        if (ovlp12>=1)
              for kk=jj+1:1:M
                  one_pos3 = zeros(1,10);
                  one_num3 = 0;
                  for oo=1:1:N
                      if (H_c_0(kk,oo)==1)
                          one_num3 = one_num3 + 1;
                          one_pos3(one_num3) = oo;
                      end
                  end
                  ovlp13 = 0;
                  ovlp13_pos = zeros(1, 12);
                  for pp=1:1:one_num1
                      for qq=1:1:one_num3
                          if (one_pos1(pp)==one_pos3(qq))
                              ovlp13 = ovlp13 + 1;
                              ovlp13_pos(ovlp13) = one_pos1(pp);
                          end
                       end
                  end
                  ovlp23 = 0;
                  ovlp23_pos = zeros(1, 12);
                  for pp=1:1:one_num2
                      for qq=1:1:one_num3
                          if (one_pos2(pp)==one_pos3(qq))
                              ovlp23 = ovlp23 + 1;
                              ovlp23_pos(ovlp23) = one_pos2(pp);
                          end
                      end
                  end
                  for pp = 1:1:ovlp12
                     for qq = 1:1:ovlp13
                         for tt=1:1:ovlp23
                         
                            if (ovlp13>=1&&ovlp23>=1&&ovlp13_pos(qq)~=ovlp23_pos(tt) && ovlp13_pos(qq)~=ovlp12_pos(pp) && ovlp23_pos(tt)~=ovlp12_pos(pp))
                               if (mod((H_c(ii,ovlp12_pos(pp))-H_c(ii,ovlp13_pos(qq))+H_c(kk,ovlp13_pos(qq))-H_c(kk,ovlp23_pos(tt))+H_c(jj,ovlp23_pos(tt))-H_c(jj,ovlp12_pos(pp))),Z)==0)
                                  g6 = g6 + 1;
                                  g6s = g6s + 1;
                                  g6_pos(g6,:) = [ii, ovlp12_pos(pp), ii, ovlp13_pos(qq), kk, ovlp13_pos(qq), kk, ...
                                  ovlp23_pos(tt), jj, ovlp23_pos(tt), jj, ovlp12_pos(pp)];
                               end
                            end
                         end
                     end
                  end
                end
            end
    end
end

                
                
                
                