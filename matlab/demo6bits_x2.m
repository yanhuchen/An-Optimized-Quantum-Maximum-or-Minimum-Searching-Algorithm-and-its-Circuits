clc
close
clear

%ȫ0̬
psi = zeros(2^6,1);
psi(1) = 1;

%������
H=1/sqrt(2)*[1,1;1,-1];
X = [0,1;1,0];
H2 = kron(H,H);
I = eye(2);
O1 = [1,0;0,0];
O2 = [0,0;0,1];

%W����������|0>�ϱ�Ϊ��̬
U_1 = kron(Ry(1.9106),eye(4));
U_2 = kron(O1,eye(4)) + kron(O2,H2);
U_3 = kron(kron(O1,H) + kron(O2,eye(2)), eye(2));
W1 = U_3*U_2*U_1;
%W1*psi%��������

U_1 = kron(Ry(1.2310),eye(4));
U_2 = kron(O1,H2) + kron(O2,eye(4));
U_3 = kron(kron(O1,eye(2)) + kron(O2,H), eye(2));
W2 = U_3*U_2*U_1;

W = kron(W1,W2);
psi = W*psi;%�ɹ��Ʊ���̬
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

%�㷨�Ļ�������
x2=43;
M2=21;
beta = asin(sqrt(M2/64));
%beta = asin(sqrt(15/36))
J = floor((pi/2-beta) / (2*beta));
phi = 2*asin(sin(pi/(4*J+6)) / sin(beta));

%oracle������д��ڵ���x0����
u1 = U3(0,0,phi);
%���48~63

C1I_3X = kron(O1,eye(2^4)) + kron(kron(O2,eye(2^3)),X);%5����
C1I_3u1 = kron(O1,eye(2^4)) + kron(kron(O2,eye(2^3)),u1);

C1_2I_3X = kron(O1,eye(2^5)) + kron(O2,C1I_3X);%6����
C1_2I_3u1 = kron(O1,eye(2^5)) + kron(O2,C1I_3u1);

%���44~47
C1IX = kron(O1,eye(4)) + kron(kron(O2,eye(2)),X);%3����
C1Iu1 = kron(O1,eye(4)) + kron(kron(O2,eye(2)),u1);

C1_2IX = kron(O1,eye(2^3)) + kron(O2,C1IX);%4����
C1_2Iu1 = kron(O1,eye(2^3)) + kron(O2,C1Iu1);

C0C1_2IX = kron(O1,C1_2IX)+ kron(O2,eye(2^4));%5����
C0C1_2Iu1= kron(O1,C1_2Iu1)+ kron(O2,eye(2^4));

C1C0C1_2IX =kron(O1,eye(2^5)) + kron(O2,C0C1_2IX);%6����
C1C0C1_2Iu1 =kron(O1,eye(2^5)) + kron(O2,C0C1_2Iu1);

%���43
C1u1 = kron(O1,eye(2)) + kron(O2,u1);%2����
C0C1u1 =  kron(O1,C1u1) + kron(O2,eye(2^2));%3����
C1C0C1u1 = kron(O1,eye(2^3)) + kron(O2,C0C1u1);%4����
C0C1C0C1u1 = kron(O1,C1C0C1u1) + kron(O2,eye(2^4));%5����
C1C0C1C0C1u1 = kron(O1,eye(2^5)) + kron(O2,C0C1C0C1u1);%6����

%���һ���޹�̬�������Գɹ����Ƿ����Ӱ��
%���17,17�����ݿ��е�һ��Ԫ��
%���29,29�������ݿ��е�Ԫ��
%���ֱ����ص�̬�ᵼ�³ɹ����½���������޹ص�̬��Ӱ��ɹ���
% C0u1 =  kron(O1,u1) + kron(O2,eye(2));%2����
% C0_2u1 = kron(O1,C0u1) + kron(O2,eye(4));%3����
% C0_3u1 = kron(O1,eye(2^3)) + kron(O2,C0_2u1);%4����
% C1C0_3u1 = kron(O1,C0_3u1) + kron(O2,eye(2^4));%5���� 
% C0C1C0_3u1 = kron(O1,C1C0_3u1) + kron(O2,eye(2^5));%6����

%���ظ����ĳ��̬�����ʲô���
%���52,52�����ݿ��е�һ��̬
%�ظ����һ��ֻ̬�е�һ�������ã��ڶ������ò���ʹ��û�б��
% C0u1 =  kron(O1,u1) + kron(O2,eye(2));%2����
% C1C0u1 = kron(O1,eye(2^2)) + kron(O2,C0u1);%3����
% C0C1C0u1 = kron(O1,C1C0u1) + kron(O2,eye(2^3));%4����
% C1C0C1C0u1 = kron(O1,eye(2^4)) + kron(O2,C0C1C0u1);%3����
% C1_2C0C1C0u1 = kron(O1,eye(2^5)) + kron(O2,C1C0C1C0u1);%3����

%���ɹ��Ƶ���̬��W���󣬼�ΪV
XRyI_4 = kron(kron(X,Ry(2.1221)),eye(2^4));

C0IRy = kron(kron(O1,eye(2)),Ry(2.2143)) + kron(O2,eye(4));%3����
IC0IRyI_2 = kron(kron(eye(2),C0IRy),eye(4));

C0X_2 = kron(O1,kron(X,X)) + kron(O2,eye(4));%3����
XC0X_2 = kron(X,C0X_2);%4bit
C0XC0X_2 = kron(O1,XC0X_2) + kron(O2,eye(2^4));%5bit
IC0XC0X_2 = kron(I,C0XC0X_2);

C1H_2 = kron(O1,eye(4)) + kron(O2,kron(H,H));%3����
XC1H_2 = kron(I,C1H_2);%4bit
C0XC1H_2 = kron(O1,XC1H_2) + kron(O2,eye(2^4));%5bit
IC0XC1H_2 = kron(I,C0XC1H_2);

H_4 = kron(kron(H,H),kron(H,H));%4bit
C1H_4 = kron(O1,eye(2^4)) + kron(O2,H_4);%5bit
IC1H_4 = kron(I,C1H_4);

V = IC1H_4*IC0XC1H_2*IC0XC0X_2*IC0IRyI_2*XRyI_4;
 
O =C1C0C1C0C1u1 * C1C0C1_2Iu1*C1C0C1_2IX*C1C0C1_2Iu1*C1C0C1_2IX * C1_2I_3u1*C1_2I_3X*C1_2I_3u1*C1_2I_3X;
%psi = O*psi;%oracle���Ӳ��Գɹ�

%������������
shift = shift_phase(6,phi);

H6 = kron(kron(kron(H,H),kron(H,H)),kron(I,I));
G = H6*shift*H6*O;%����W���ڳ���W����
for i=1:J+1
    psi = G*psi;
end
abs(psi(63))^2*15