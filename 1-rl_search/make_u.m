clc

fp = fopen('u_tab.dat', 'w');
for ii=1:1:10001
    fprintf(fp, ' %d', u_tab(ii));
end
fclose(fp);