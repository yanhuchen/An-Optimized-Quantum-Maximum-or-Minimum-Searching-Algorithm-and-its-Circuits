clc
close
clear

psi = zeros(2^3,1);
psi(1) = 1;

H=1/sqrt(2)*[1,1;1,-1];
X = [0,1;1,0];
H2 = kron(H,H);
O1 = [1,0;0,0];
O2 = [0,0;0,1];

% 已知有初态包含了7个数据，最后一位为0，其余为1/sqrt(7)，
% 若第一次随机选择的x0=4，则需要标记4,5,6,7,。但是因为一共只有前7个有值
% 所以标记的时候只能标记4,5,6三个态。sin(β)=sqrt(3/7)
% 在得到这个结论之前，我们已经尝试过 sin(β)=sqrt(3/8)和sin(β)=sqrt(4/8)两种情况
% 最直接的是将大于等于x0的所有态全部标记出来，但是这样会有一定误差

%phi=pi/2;
phi = 1.7382
%Oracle作用是令状态4,5,6,7添加一个相位
CIX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);
u1 = U3(0,0,phi);
CIu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);
C1C0u1 =  kron(O1,eye(4)) + kron(O2, kron(O1,u1)+kron(O2,eye(2)));
O = CIu1*CIX*CIu1*CIX;%标记4,5,6，7
%O = C1C0u1*CIX*CIu1*CIX;%标记4,5,6
%条件相移算子
C0X = kron(O1,X) + kron(O2,eye(2));
C0_2X = kron(O1,C0X) + kron(O2,eye(4));
C0u1 = kron(O1,u1) + kron(O2,eye(2));
C0_2u1 = kron(O1,C0u1) + kron(O2,eye(4));
shift = C0_2X*C0_2u1*C0_2X;

%W算子
U1 = kron(Ry(1.427444876),eye(4));
U2 = kron(O1,H2) + kron(O2,eye(4));
U3 = kron(kron(O1,eye(2)) +  kron(O2,Ry(1.23095942)),eye(2));
U4 = kron(O1,eye(4)) + kron(O2 ,kron(O1,H) + kron(O2,eye(2)))
%U5 = kron(O1,eye(4)) + kron(O2,kron(O1,eye(2)) + kron(O2,X));
W =  U4*U3*U2*U1;

state0 = W*psi
state1 = O*state0
state2 = inv(W)*state1
state3 = shift*state2
state4 = W*state3