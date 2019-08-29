clc
close
clear

%基础数据
N=25;
n=ceil(log2(N));
x0=13;
M0 = N-x0;
beta = asin(sqrt(M0/2^n));
J = floor((pi/2-beta)/(2*beta));
phi = 2*asin((pi/(4*J+6))/sin(beta));

%基本门
I = eye(2);
X = [0,1;1,0];
H = 1/sqrt(2)*[1,1;1,-1];
H4 = kron(kron(H,H),kron(H,H));
H5  = kron(H4,H);
%用于构建受控门的算子
O1 = [1,0;0,0];
O2 = [0,0;0,1];
%u1门
u1 = U3(0,0,phi);

%各种受控门
C1I_3X = kron(O1,eye(2^4))+kron(kron(O2,eye(2^3)),X);
C1I_3u1 = kron(O1,eye(2^4))+kron(kron(O2,eye(2^3)),u1);

C1X = kron(O1,eye(2)) + kron(O2,X);
C1_2X = kron(O1,eye(4)) + kron(O2,C1X);
C1_3X = kron(O1,eye(8)) + kron(O2,C1_2X);
C0C1_3X = kron(O1,C1_3X) + kron(O2,eye(16));

C1u1 = kron(O1,eye(2)) + kron(O2,u1);
C1_2u1 = kron(O1,eye(4)) + kron(O2,C1u1);
C1_3u1 = kron(O1,eye(8)) + kron(O2,C1_2u1);
C0C1_3u1 = kron(O1,C1_3u1) + kron(O2,eye(16));

C1Iu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);
C1_2Iu1 = kron(O1,eye(8)) + kron(O2,C1Iu1);
C0C1_2Iu1 = kron(O1,C1_2Iu1) + kron(O2,eye(16));

C0X = kron(O1,X) + kron(O2,eye(2));
C0_2X = kron(O1,C0X) + kron(O2,eye(4));
C0_3X = kron(O1,C0_2X) + kron(O2,eye(8));
C0_4X = kron(O1,C0_3X) + kron(O2,eye(16));

C0u1 = kron(O1,u1) + kron(O2,eye(2));
C0_2u1 = kron(O1,C0u1) + kron(O2,eye(4));
C0_3u1 = kron(O1,C0_2u1) + kron(O2,eye(8));
C0_4u1 = kron(O1,C0_3u1) + kron(O2,eye(16));

%设计初态|psi>
psi = zeros(32,1);
psi(1:25)=1/5;

%Oracle设计
U1 = C1I_3u1*C1I_3X*C1I_3u1*C1I_3X;%标记16—31
U2 = C0C1_3X*C0C1_3u1*C0C1_3X;%标记14
U3 = C0C1_2Iu1;%标记13、15
O= U3*U2*U1;

%条件相移设计
shift =C0_4X+C0_4u1*C0_4X; 

for i=1:J+1
    psi = H5*shift*H5*O*psi;
end