clc
close
clear

%全0态
psi = zeros(2^6,1);
psi(1) = 1;

%基本门
H=1/sqrt(2)*[1,1;1,-1];
X = [0,1;1,0];
H2 = kron(H,H);
I = eye(2);
O1 = [1,0;0,0];
O2 = [0,0;0,1];

%W算子作用在|0>上变为初态
U_1 = kron(Ry(1.9106),eye(4));
U_2 = kron(O1,eye(4)) + kron(O2,H2);
U_3 = kron(kron(O1,H) + kron(O2,eye(2)), eye(2));
W1 = U_3*U_2*U_1;
%W1*psi%测试无误

U_1 = kron(Ry(1.2310),eye(4));
U_2 = kron(O1,H2) + kron(O2,eye(4));
U_3 = kron(kron(O1,eye(2)) + kron(O2,H), eye(2));
W2 = U_3*U_2*U_1;

W = kron(W1,W2);
psi = W*psi;%成功制备初态
% k=1;
% for i=1:64
%     if psi(i) ~=0
%         R(k) = i-1;
%         k=k+1;
%     end
% end
% index = randperm(36);
% for i =1:36
%     R1(i) = R(index(i));
% end

%算法的基本参数
x2=43;
M2=21;
beta = asin(sqrt(M2/64));
%beta = asin(sqrt(15/36))
J = floor((pi/2-beta) / (2*beta));
phi = 2*asin(sin(pi/(4*J+6)) / sin(beta));

%oracle标记所有大于等于x0的数
u1 = U3(0,0,phi);
%标记48~63

C1I_3X = kron(O1,eye(2^4)) + kron(kron(O2,eye(2^3)),X);%5比特
C1I_3u1 = kron(O1,eye(2^4)) + kron(kron(O2,eye(2^3)),u1);

C1_2I_3X = kron(O1,eye(2^5)) + kron(O2,C1I_3X);%6比特
C1_2I_3u1 = kron(O1,eye(2^5)) + kron(O2,C1I_3u1);

%标记44~47
C1IX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);%3比特
C1Iu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);

C1_2IX = kron(O1,eye(2^3)) + kron(O2,C1IX);%4比特
C1_2Iu1 = kron(O1,eye(2^3)) + kron(O2,C1Iu1);

C0C1_2IX = kron(O1,C1_2IX)+ kron(O2,eye(2^4));%5比特
C0C1_2Iu1= kron(O1,C1_2Iu1)+ kron(O2,eye(2^4));

C1C0C1_2IX =kron(O1,eye(2^5)) + kron(O2,C0C1_2IX);%6比特
C1C0C1_2Iu1 =kron(O1,eye(2^5)) + kron(O2,C0C1_2Iu1);

%标记43
C1u1 = kron(O1,eye(2)) + kron(O2,u1);%2比特
C0C1u1 =  kron(O1,C1u1) + kron(O2,eye(2^2));%3比特
C1C0C1u1 = kron(O1,eye(2^3)) + kron(O2,C0C1u1);%4比特
C0C1C0C1u1 = kron(O1,C1C0C1u1) + kron(O2,eye(2^4));%5比特
C1C0C1C0C1u1 = kron(O1,eye(2^5)) + kron(O2,C0C1C0C1u1);%6比特

%标记一个无关态，看看对成功率是否存在影响
%标记17,17是数据库中的一个元素
%标记29,29不是数据库中的元素
%发现标记相关的态会导致成功率下降，而标记无关的态则不影响成功率
% C0u1 =  kron(O1,u1) + kron(O2,eye(2));%2比特
% C0_2u1 = kron(O1,C0u1) + kron(O2,eye(4));%3比特
% C0_3u1 = kron(O1,eye(2^3)) + kron(O2,C0_2u1);%4比特
% C1C0_3u1 = kron(O1,C0_3u1) + kron(O2,eye(2^4));%5比特 
% C0C1C0_3u1 = kron(O1,C1C0_3u1) + kron(O2,eye(2^5));%6比特

%若重复标记某个态会出现什么结果
%标记52,52是数据库中的一个态
%重复标记一个态只有第一次起作用，第二次作用不会使它没有标记
% C0u1 =  kron(O1,u1) + kron(O2,eye(2));%2比特
% C1C0u1 = kron(O1,eye(2^2)) + kron(O2,C0u1);%3比特
% C0C1C0u1 = kron(O1,C1C0u1) + kron(O2,eye(2^3));%4比特
% C1C0C1C0u1 = kron(O1,eye(2^4)) + kron(O2,C0C1C0u1);%3比特
% C1_2C0C1C0u1 = kron(O1,eye(2^5)) + kron(O2,C1C0C1C0u1);%3比特

%生成估计叠加态的W矩阵，记为V
XRyI_4 = kron(kron(X,Ry(2.1221)),eye(2^4));

C0IRy = kron(kron(O1,eye(2)),Ry(2.2143)) + kron(O2,eye(4));%3比特
IC0IRyI_2 = kron(kron(eye(2),C0IRy),eye(4));

C0X_2 = kron(O1,kron(X,X)) + kron(O2,eye(4));%3比特
XC0X_2 = kron(X,C0X_2);%4bit
C0XC0X_2 = kron(O1,XC0X_2) + kron(O2,eye(2^4));%5bit
IC0XC0X_2 = kron(I,C0XC0X_2);

C1H_2 = kron(O1,eye(4)) + kron(O2,kron(H,H));%3比特
XC1H_2 = kron(I,C1H_2);%4bit
C0XC1H_2 = kron(O1,XC1H_2) + kron(O2,eye(2^4));%5bit
IC0XC1H_2 = kron(I,C0XC1H_2);

H_4 = kron(kron(H,H),kron(H,H));%4bit
C1H_4 = kron(O1,eye(2^4)) + kron(O2,H_4);%5bit
IC1H_4 = kron(I,C1H_4);

V = IC1H_4*IC0XC1H_2*IC0XC0X_2*IC0IRyI_2*XRyI_4;
 
O =C1C0C1C0C1u1 * C1C0C1_2Iu1*C1C0C1_2IX*C1C0C1_2Iu1*C1C0C1_2IX * C1_2I_3u1*C1_2I_3X*C1_2I_3u1*C1_2I_3X;
%psi = O*psi;%oracle算子测试成功

%条件相移算子
shift = shift_phase(6,phi);

H6 = kron(kron(kron(H,H),kron(H,H)),kron(I,I));
G = H6*shift*H6*O;%除以W等于乘以W的逆
for i=1:J+1
    psi = G*psi;
end
abs(psi(63))^2*15