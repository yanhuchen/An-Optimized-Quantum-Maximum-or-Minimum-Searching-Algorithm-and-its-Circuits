clc
close
clear 

N=1024;
mu=N/2;%均值
sigma=N/4;%标准差
x = mu-2*sigma:mu+2*sigma;
x = x';
%高斯函数
fx = 1/(sqrt(2*pi)*sigma) * exp(-(x-mu).^2 / (2*sigma^2) );
%plot(x,fx)
t=sum(fx)

rk = zeros(N+1,1);
for j=1:N+1
    for i= 1:N+1
        if sum(fx(1:i))/t > (j-1)/(N+1)
            rk(j) = i;
            break;
        end
    end
end
x = x/N
rk = rk/N
p = polyfit(x,rk,4);

f = polyval(p,x);
plot(x,rk,'b-',x,f,'r-')

