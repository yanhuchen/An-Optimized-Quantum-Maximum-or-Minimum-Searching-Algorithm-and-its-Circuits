clc
close 
clear
[x,y]=meshgrid(0:1:4,0:1:13);%用户的位置
% sinb1=sqrt(x);
% sinb2=sqrt(y);
% cosb1=sqrt(1-sinb1.*sinb1);
%  b2=asin(sinb2);
% 
% for k=1:1:length(x)
%     for l=1:1:length(y)
% 
%         J=floor(((pi/2)-b2(l,1))/(2*b2(l,1)));
%         %disp(J);
%         F=2*asin(sin(pi/(4*J+6))/sinb2(l,1));
%         %disp(F);
%         E=exp(1i*F);
%       % disp(E);
%         A=[-E*(1+(E-1)*sinb1(1,k)*sinb1(1,k)),-(E-1)*sinb1(1,k)*cosb1(1,k);-E*(E-1)*sinb1(1,k)*cosb1(1,k),-E+(E-1)*sinb1(1,k)*sinb1(1,k)];
%         %disp(A);
%         B=[sinb1(1,k);cosb1(1,k)];
%         %disp(B);
%         D=A^(J+1);
%         %disp(D);
%         X=D*B;
%         %disp(X);
%          z(l,k)=abs(X(1,1))*abs(X(1,1));
% %z=abs(X(1,1))*abs(X(1,1));
%     end
% end
% % [m,n]=size(X);
% % 
% % disp(z);
% sum1=0;
% for q=1:1:length(x)
%     for t=1:1:q
%        sum1=sum1+z(q,t);       
%     end
% end
% disp(sum1);
% sum2=0;
% for m=1:1:length(x)
%     for n=m:1:length(x)
%        sum2=sum2+z(m,n);       
%     end
% end
% disp(sum2);

figure (1)
z=[0.044 0.052 0.279 0.629 0;
   0.045 0.041 0.479 0.436 0;
   0.066 0.222 0.340 0.372 0;
   0.089 0.160 0.330 0.422 0;
   0.113 0.097 0.134 0.656 0;
   0.107 0.127 0.139 0.627 0;
   0.075 0.050 0.702 0.173 0;
   0.060 0.100 0.362 0.478 0;
   0.085 0.169 0.241 0.504 0;
   0.095 0.124 0.301 0.480 0;
   0.103 0.242 0.296 0.357 0;
   0.144 0.211 0.401 0.244 0;
   0,0,0,0,0;
   0,0,0,0,0]
surf(0:1:4,0:1:13,z);%用surf画出来的图叫surf三维着色表面图，x,y通常是通过调用meshgrid函数生成的数据网络
view(0,90)
xlabel('M1/N (真实)');ylabel('M2/N (假设)');zlabel('P(成功率)');   
