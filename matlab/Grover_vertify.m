clc
close
clear

H=1/sqrt(2)*[1,1;1,-1];
I=eye(2);
W = kron(kron(H,H),I);
I0 = -eye(4);
I0(1) = 1;
I0 = kron(I0,I);
It = eye(8);
It(7,7) = 0;
It(7,8) = 1;
It(8,8) = 0;
It(8,7) = 1

G = -W*I0*W*It
state0 = 1/2/sqrt(2)*[1;-1;1;-1;1;-1;1;-1];
G*state0
N=8;
beta = asin(sqrt(1/N))
J = floor((pi/2-beta)/(2*beta))
phi = 2*asin(sin(pi/(4*J+6))/sqrt(1/N))