clc
close
clear

psi = zeros(4,1);
psi(1) = 1;

H=1/sqrt(2)*[1,1;1,-1];
X = [0,1;1,0];
H2 = kron(H,H);
O1 = [1,0;0,0];
O2 = [0,0;0,1];
phi=pi/2;
%Oracle��������״̬2,3���һ����λ