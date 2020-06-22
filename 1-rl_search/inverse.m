%clear all
clc

%v=1:-0.0001:0;
u1 = zeros(1,10001);
temp = zeros(1,10001);
tp = 1;
for ii=1:1:10001
    tp =1;
    for jj=tp:1:1000001
        if (y_o(jj) < (1-(ii-1)*0.0001))
            u1(ii) = (jj-2)*0.0001;
            temp(ii) = jj-1;
            tp=jj;
            break;
        end
    end
end
u1(1) = 0;
u1(10001) = 2*u1(10000)-u1(9999);

u_tab = zeros(1,10001);
for ii=1:1:10001
    u_tab(ii) = u1(10002-ii);
end

err_num = 0;
for ii=2:1:10001
    if (u_tab(ii)==u_tab(ii-1))
        err_num = err_num + 1;
        err(err_num) = ii;
    end
end

y=zeros(1, 35001);
for ii =1:1:35001
    y(ii) = y_o((ii-1)*10 + 1);
end

v = 0:0.0001:1;
plot(v, u_tab);
hold on
x=0:0.001:35;
plot(x, y);
grid on
