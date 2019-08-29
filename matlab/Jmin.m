clc
close
clear

 N=2;%����1��
% N=200;%����10��
% N=16500;%����100��
% N=1625000;%����1000��
M=1;
psi = [sqrt(M/N);sqrt((N-M)/N)];%׼ȷ��̬ʸ��
iter = 2000;
epi = zeros(iter,1);
J = zeros(iter,1);
%ʹ�ô�������̬ʸ�������beta��J��phi��G��ȫ�����Ǵ������Ĺ���ֵ
%����Grover������G���ƾ������׼ȷ��̬ʸ��������ģ��ʵ������У���֪����;�����̬����
for i =1:iter
   epi(i) = (i-1)/iter*0.3;
   psi_pie = [sqrt(M/N + epi(i));
              sqrt((N-M)/N - epi(i))];%��������̬ʸ��
   beta = asin(psi_pie(1));
   J(i) = floor((pi/2-beta) / (2*beta));
   while sin(pi/(4*J+6)) > sin(beta)
       J(i) = J(i)+1;
   end
   %û�������ӳɹ��ʵ�����
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