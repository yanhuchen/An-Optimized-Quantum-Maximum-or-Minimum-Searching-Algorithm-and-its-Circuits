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
%c0c0x�ţ���ʾ��00ʱ�ܿر��ط�ת
c0c0x = kron(O1,c0x) +  kron(O2,eye(4));
c0c0z = kron(O1,c0z) +  kron(O2,eye(4));
%c1_3x��ʾ��111ʱ���ܿر��ط�ת
c1_3x =  kron(O1,eye(8)) + kron(O2,c1c1x);


%һ��3bit��Grover����
%oracle����
O = c1_3x;%���̬Ϊ7

%Hn
H3I = kron(kron(kron(H,H),H),I);
%������������
shift = kron(c0c0x*c0c0z*c0c0x,I);
%ע���4���������һ���ʱ���˳��
G = H3I * shift * H3I * O;

state_0 = eye(2^4,1)%ȫ0̬
%��ʼ�Ĵ���֮���Ϊ��̬����Grvoer�㷨��ֻ���ʼ����һ��
start = kron(kron(kron(H,H),H),H*X);
%������ɺ����β������ֻ��������һ��
over = kron(eye(8),H);
state_end = over*G*start*state_0

