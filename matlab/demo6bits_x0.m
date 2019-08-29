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
x0=19;
M0=45;
beta = asin(sqrt(M0/64));
J = floor((pi/2-beta) / (2*beta));
phi = 2*asin(sin(pi/(4*J+6)) / sin(beta));

%oracle标记所有大于等于x0的数
u1 = U3(0,0,phi);
%标记32~63
CI_4X = kron(O1,eye(2^5)) + kron(kron(O2,eye(2^4)),X);%6比特
CI_4u1 = kron(O1,eye(2^5)) + kron(kron(O2,eye(2^4)),u1);

C1IX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);%3比特
C1Iu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);

IC1IX = kron(eye(2),C1IX);%4比特
IC1Iu1 = kron(eye(2),C1Iu1);

C1IC1IX = kron(O1,eye(2^4))+ kron(O2,IC1IX);%5比特
C1IC1Iu1 = kron(O1,eye(2^4))+ kron(O2,IC1Iu1);
%标记20~31
C0C1IC1IX = kron(O1,C1IC1IX)+ kron(O2,eye(2^5));%6比特
C0C1IC1Iu1 = kron(O1,C1IC1Iu1)+ kron(O2,eye(2^5));
%标记19
C1u1 = kron(O1,eye(2)) + kron(O2,u1);%2比特
C0C1u1 = kron(O1,C1u1) + kron(O2,eye(2^2));%3比特
C0_2C1u1 = kron(O1,C0C1u1) + kron(O2,eye(2^3));%4比特
C1C0_2C1u1 = kron(O1,eye(2^4)) + kron(O2,C0_2C1u1);%5比特
C0C1C0_2C1u1 = kron(O1,C1C0_2C1u1) + kron(O2,eye(2^5));%6比特

O = C0C1C0_2C1u1 * C0C1IC1Iu1*C0C1IC1IX*C0C1IC1Iu1*C0C1IC1IX * CI_4u1*CI_4X*CI_4u1*CI_4X;
%psi = O*psi;%oracle算子测试成功

%条件相移算子
shift = shift_phase(6,phi);
G = W*shift/W*O;
for i=1:J+1
    psi = G*psi;
end

