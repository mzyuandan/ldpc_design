function flag = girth_cmp(girth, set)

N = length(girth);

flag = 1;
for ii=1:1:N
    if (girth(ii) > set(ii))
        flag = 0;
        break;
    end
end

end
