function [U_t] = shfit_phase(n,phi)
X = [0,1;1,0];
O1 = [1,0;0,0];
O2 = [0,0;0,1];
u1 = U3(0,0,phi);
U=X;
U_1=u1;
for i=2:n
    U = kron(O1,U) + kron(O2,eye(2^(i-1)));
    U_1 = kron(O1,U_1) + kron(O2,eye(2^(i-1)));
end
U_t = U*U_1*U;
end