clc
close
clear

psi = [1;0;0;0;0;0;0;0];
H=1/sqrt(2)*[1,1;
             1,-1];
X = [0,1;0,1];
O1 = [1,0;0,0];
O2 = [0,0;0,1];
I = eye(2);
RyII = kron(Ry(1.2310),eye(4));
C0HI = kron(kron(O1,H) + kron(O2,I),I);
C0IH = kron(kron(O1,I),H) + kron(O2,eye(4));
C1HI = kron(O1,eye(4)) + kron(kron(O2,H),I);
%0,1,2,3,4,6¹²6¸öÌ¬
W = C1HI*C0IH*C0HI*RyII*psi

O