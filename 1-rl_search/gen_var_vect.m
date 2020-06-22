function [ var_vect, edg_num ] = gen_var_vect( mb, nb, dv )
%gen_var_vect Summary of this function goes here
%   Detailed explanation goes here

var_vect = zeros(1, dv-1);

nv_rang = nb+1;
% while var_vect(1) < mb-1 || var_vect(1)>nb-1
%     var_vect(1) = randint(1,1,nv_rang);
% end
var_vect(1) = mb-1;

nv_rang = nv_rang - var_vect(1);
while var_vect(2) < 1
    var_vect(2) = randint(1,1,nv_rang);
end

for jj=4:1:dv-1
    nv_rang = nv_rang - var_vect(jj-2);
    var_vect(jj-1) = randint(1,1,nv_rang);
end

nv_rang = nv_rang - var_vect(dv-2);
var_vect(dv-1) = nv_rang-1;

edg_num = 0;
for jj=2:1:dv
    edg_num = edg_num + jj*var_vect(jj-1);
end

end