function [ chk_vect, err ] = gen_chk_vect(mb, edg_num, dc)
%gen_var_vect Summary of this function goes here
%   Detailed explanation goes here

test_num = 0;
err = 0;
while (1)
chk_vect = zeros(1, 3);
nc_rang = mb+1;
chk_vect(1) = randint(1,1,nc_rang);
edg_used = (dc-2) * chk_vect(1); 
while edg_used > edg_num
    chk_vect(1) = randint(1,1,nc_rang);
    edg_used = (dc-2) * chk_vect(1); 
end

nc_rang = nc_rang - chk_vect(1);
chk_vect(2) = randint(1,1,nc_rang);
edg_used_r = edg_used;
edg_used = edg_used + (dc-1) * chk_vect(2); 
while edg_used > edg_num
    chk_vect(2) = randint(1,1,nc_rang);
    edg_used = edg_used_r + (dc-1) * chk_vect(2); 
end

nc_rang = nc_rang - chk_vect(2);
chk_vect(3) = nc_rang - 1;
edg_used = edg_used + (dc) * chk_vect(3); 

if (edg_used==edg_num)
    break;
end
test_num = test_num + 1;
if (test_num >1000)
   err = 1;
   break;
end

end

end