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

% ��֪�г�̬������7�����ݣ����һλΪ0������Ϊ1/sqrt(7)��
% ����һ�����ѡ���x0=4������Ҫ���4,5,6,7,��������Ϊһ��ֻ��ǰ7����ֵ
% ���Ա�ǵ�ʱ��ֻ�ܱ��4,5,6����̬��sin(��)=sqrt(3/7)
% �ڵõ��������֮ǰ�������Ѿ����Թ� sin(��)=sqrt(3/8)��sin(��)=sqrt(4/8)�������
% ��ֱ�ӵ��ǽ����ڵ���x0������̬ȫ����ǳ�����������������һ�����

%phi=pi/2;
phi = 1.7382
%Oracle��������״̬4,5,6,7���һ����λ
CIX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);
u1 = U3(0,0,phi);
CIu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);
C1C0u1 =  kron(O1,eye(4)) + kron(O2, kron(O1,u1)+kron(O2,eye(2)));
O = CIu1*CIX*CIu1*CIX;%���4,5,6��7
%O = C1C0u1*CIX*CIu1*CIX;%���4,5,6
%������������
C0X = kron(O1,X) + kron(O2,eye(2));
C0_2X = kron(O1,C0X) + kron(O2,eye(4));
C0u1 = kron(O1,u1) + kron(O2,eye(2));
C0_2u1 = kron(O1,C0u1) + kron(O2,eye(4));
shift = C0_2X*C0_2u1*C0_2X;

%W����
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