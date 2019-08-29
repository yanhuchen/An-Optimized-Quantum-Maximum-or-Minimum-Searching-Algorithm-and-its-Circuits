clc
close 
clear

%3比特的非线性叠加态测试
psi = zeros(2^3,1);
psi(1) = 1;
H=1/sqrt(2)*[1,1;1,-1];
X = [0,1;1,0];
H2 = kron(H,H);
O1 = [1,0;0,0];
O2 = [0,0;0,1];

U1 = kron(Ry(1.4274),eye(4));
U2 = kron(O1,H2) + kron(O2,eye(4))
U3 = kron(kron(O1,eye(2)) +  kron(O2,Ry(1.2310)),eye(2));
U4 = kron(O1,eye(4)) + kron(O2 ,kron(O1,H) + kron(O2,eye(2)))
W =  U4*U3*U2*U1
