function [U] = U3(theta,phi,lamda)
U=[cos(theta/2),  -exp(1i*lamda)*sin(theta/2);
   exp(1i*phi)*sin(theta/2),  exp(1i*(phi+lamda))*cos(theta/2)];


end

