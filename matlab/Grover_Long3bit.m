clc
close 
clear

phi = 2.1269
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
%c0c0x�ţ���ʾ��00ʱ�ܿر��ط�ת
c0c0x = kron(O1,c0x) +  kron(O2,eye(4));
c0c0z = kron(O1,c0z) +  kron(O2,eye(4));
%�ܿ�u3��
cu3 = kron(O1,I)+kron(O2,U3(pi,2.1269,-1.0147));
ccu3 = kron(O1,eye(4))+kron(O2,cu3);
c3u3 = kron(O1,eye(8))+kron(O2,ccu3);
%�ܿ�u1��
u1 = U3(0,0,2.1269);
c0u1 = kron(O1,u1) +  kron(O2,eye(2));
c0c0u1 = kron(O1,c0u1) +  kron(O2,eye(4));
%һ��3bit��Grover����
%oracle����
O = c3u3;%���̬Ϊ7

%Hn
H3I = kron(kron(kron(H,H),H),I);
%������������
shift = kron(c0c0x*c0c0u1*c0c0x,I);
%ע���4���������һ���ʱ���˳��
G = H3I * shift * H3I * O;

state_0 = eye(2^4,1);%ȫ0̬
%��ʼ�Ĵ���֮���Ϊ��̬����Grvoer�㷨��ֻ���ʼ����һ��
start = kron(kron(kron(H,H),H),H*X);
%������ɺ����β������ֻ��������һ��
over = kron(eye(8),H);
state_1 = start*state_0
state_2 = O*start*state_0
state_3 = H3I*O*start*state_0
state_4 = shift*H3I*O*start*state_0
state_5 = H3I*shift*H3I*O*start*state_0
state_6 = over*H3I*shift*H3I*O*start*state_0
state_end = over*G*G*G*start*state_0