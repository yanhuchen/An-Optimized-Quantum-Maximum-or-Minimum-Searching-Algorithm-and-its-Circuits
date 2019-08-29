clc
close 
clear

%基本门
I = eye(2);
X = [0,1;1,0];
H = 1/sqrt(2)*[1,1;1,-1];
Z = [1,0;0,-1];
%用于构建受控门的算子
O1 = [1,0;0,0];
O2 = [0,0;0,1];

%各种受控门
c0x = kron(O1,X) + kron(O2,I);
c1x = kron(O1,I) + kron(O2,X);
c0z = kron(O1,Z) + kron(O2,I);
%c1c1x门，表示在11时受控比特翻转
c1c1x = kron(O1,eye(4)) + kron(O2,c1x);

%一次2bit的Grover迭代
%oracle算子
O = c1c1x;%标记态为3
%Hn
H2I = kron(kron(H,H),I);
%条件相移算子
shift = kron(c0x*c0z*c0x,I);
%注意把4个步骤合在一起的时候的顺序
G = H2I * shift * H2I * O;

state_0 = eye(2^3,1)%全0态
%初始的处理，之后变为初态输入Grvoer算法，只在最开始运行一次
start = kron(kron(H,H),H*X);
%迭代完成后的收尾工作，只在最后进行一次
over = kron(eye(4),H);
state_1 = start*state_0
state_2 = O*start*state_0
state_3 = H2I*O*start*state_0
state_4 = shift*H2I*O*start*state_0
state_5 = H2I*shift*H2I*O*start*state_0
state_5_pie = G*start*state_0
state_end = over*G*start*state_0