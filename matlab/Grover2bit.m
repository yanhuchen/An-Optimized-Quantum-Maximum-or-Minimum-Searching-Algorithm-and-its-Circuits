clc
close 
clear

%������
I = eye(2);
X = [0,1;1,0];
H = 1/sqrt(2)*[1,1;1,-1];
Z = [1,0;0,-1];
%���ڹ����ܿ��ŵ�����
O1 = [1,0;0,0];
O2 = [0,0;0,1];

%�����ܿ���
c0x = kron(O1,X) + kron(O2,I);
c1x = kron(O1,I) + kron(O2,X);
c0z = kron(O1,Z) + kron(O2,I);
%c1c1x�ţ���ʾ��11ʱ�ܿر��ط�ת
c1c1x = kron(O1,eye(4)) + kron(O2,c1x);

%һ��2bit��Grover����
%oracle����
O = c1c1x;%���̬Ϊ3
%Hn
H2I = kron(kron(H,H),I);
%������������
shift = kron(c0x*c0z*c0x,I);
%ע���4���������һ���ʱ���˳��
G = H2I * shift * H2I * O;

state_0 = eye(2^3,1)%ȫ0̬
%��ʼ�Ĵ���֮���Ϊ��̬����Grvoer�㷨��ֻ���ʼ����һ��
start = kron(kron(H,H),H*X);
%������ɺ����β������ֻ��������һ��
over = kron(eye(4),H);
state_1 = start*state_0
state_2 = O*start*state_0
state_3 = H2I*O*start*state_0
state_4 = shift*H2I*O*start*state_0
state_5 = H2I*shift*H2I*O*start*state_0
state_5_pie = G*start*state_0
state_end = over*G*start*state_0