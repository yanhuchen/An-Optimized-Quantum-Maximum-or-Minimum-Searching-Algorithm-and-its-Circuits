clc
close
clear

 N=2;%迭代1次
% N=200;%迭代10次
% N=16500;%迭代100次
% N=1625000;%迭代1000次
M=1;
psi = [sqrt(M/N);sqrt((N-M)/N)];%准确的态矢量
iter = 2000;
epi = zeros(iter,1);
J = zeros(iter,1);
%使用带有误差的态矢量算出的beta，J，phi，G，全部都是带有误差的估计值
%但是Grover迭代的G估计矩阵乘以准确的态矢量，用以模拟实际情况中，不知排序和具体标记态数。
for i =1:iter
   epi(i) = (i-1)/iter*0.3;
   psi_pie = [sqrt(M/N + epi(i));
              sqrt((N-M)/N - epi(i))];%带有误差的态矢量
   beta = asin(psi_pie(1));
   J(i) = floor((pi/2-beta) / (2*beta));
   while sin(pi/(4*J+6)) > sin(beta)
       J(i) = J(i)+1;
   end
   %没有起到增加成功率的作用
   %J(i) = J(i)+2;
   phi= 2*asin(sin(pi/(4*J(i)+6)) / sin(beta));
   
   G = [-exp(1i*phi) * (1+(exp(1i*phi)-1) * sin(beta)*sin(beta)),  -(exp(1i*phi)-1) *sin(beta)*cos(beta);
     -exp(1i*phi) * (exp(1i*phi)-1) * sin(beta)*cos(beta),   -exp(1i*phi) + (exp(1i*phi)-1) * sin(beta)*sin(beta)];
   a =  G^(J(i)+1)*psi;
   solu(i) = a(2);
   solu(i) = abs(solu(i))^2;
end
solu = solu';
plot(1:iter,solu)