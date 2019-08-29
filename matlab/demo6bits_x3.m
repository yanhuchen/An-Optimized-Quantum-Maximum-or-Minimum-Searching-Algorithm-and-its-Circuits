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
x2=58;
M3=6;
beta = asin(sqrt(M3/64));
J = floor((pi/2-beta) / (2*beta));
phi = 2*asin(sin(pi/(4*J+6)) / sin(beta));

%oracle标记所有大于等于x0的数
u1 = U3(0,0,phi);

%标记60~63
C1IX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);%3比特
C1Iu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);

C1_2IX = kron(O1,eye(2^3)) + kron(O2,C1IX);%4比特
C1_2Iu1 = kron(O1,eye(2^3)) + kron(O2,C1Iu1);

C1_3IX = kron(O1,eye(2^4)) + kron(O2,C1_2IX);%5比特
C1_3Iu1 = kron(O1,eye(2^4)) + kron(O2,C1_2Iu1);

C1_4IX = kron(O1,eye(2^5)) + kron(O2,C1_3IX);%6比特
C1_4Iu1 = kron(O1,eye(2^5)) + kron(O2,C1_3Iu1);

%标记58~59
C1X = kron(O1,eye(2)) + kron(O2,X);%2比特
C1u1 = kron(O1,eye(2)) + kron(O2,u1);

C0C1IX = kron(O1,C1X) + kron(O2,eye(4));%3比特
C0C1Iu1 = kron(O1,C1u1) + kron(O2,eye(4));

C1C0C1IX = kron(O1,eye(2^3)) + kron(O2,C0C1IX);%4比特
C1C0C1Iu1 = kron(O1,eye(2^3)) + kron(O2,C0C1Iu1);

C1_2C0C1IX = kron(O1,eye(2^4)) + kron(O2,C1C0C1IX);%5比特
C1_2C0C1Iu1 = kron(O1,eye(2^4)) + kron(O2,C1C0C1Iu1);

C1_3C0C1IX = kron(O1,eye(2^5)) + kron(O2,C1_2C0C1IX);%6比特
C1_3C0C1Iu1 = kron(O1,eye(2^5)) + kron(O2,C1_2C0C1Iu1);

O = C1_3C0C1Iu1 * C1_3C0C1IX*C1_3C0C1Iu1*C1_3C0C1IX * C1_4Iu1*C1_4IX*C1_4Iu1*C1_4IX; 
%psi = O*psi;%oracle算子测试成功

%条件相移算子
shift = shift_phase(6,phi);
G = W*shift*inv(W)*O;
for i=1:J+1
    psi = G*psi;
end