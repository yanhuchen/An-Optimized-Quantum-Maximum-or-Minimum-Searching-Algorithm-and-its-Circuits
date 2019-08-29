clc
close
clear


e = 1000;
psi = [sqrt(1/4);sqrt(3/4)];
psi_pie = zeros(2,e);
beta = asin(psi(1));
J = floor((pi/2-beta) / (2*beta));
beta_pie = zeros(e,1);
phi_pie = zeros(e,1);
solu = zeros(e,1);
for i =1:e
    
    beta_pie(i) = beta + i*0.1*pi/e;
    phi_pie(i)= 2*asin(sin(pi/(4*J+6)) / sin(beta_pie(i)));
    G = [-exp(1i*phi_pie(i)) * (1+(exp(1i*phi_pie(i))-1) * sin(beta_pie(i))*sin(beta_pie(i))),  -(exp(1i*phi_pie(i))-1) *sin(beta_pie(i))*cos(beta_pie(i));
         -exp(1i*phi_pie(i)) * (exp(1i*phi_pie(i))-1) * sin(beta_pie(i))*cos(beta_pie(i)),   -exp(1i*phi_pie(i)) + (exp(1i*phi_pie(i))-1) * sin(beta_pie(i))*sin(beta_pie(i))];
    a = G^(J+1)*psi;
    solu(i) = a(2);
    solu(i) = abs(solu(i))^2;
end

phi= 2*asin(sin(pi/(4*J+6)) / sin(beta));
G = [-exp(1i*phi) * (1+(exp(1i*phi)-1) * sin(beta)*sin(beta)),  -(exp(1i*phi)-1) *sin(beta)*cos(beta);
     -exp(1i*phi) * (exp(1i*phi)-1) * sin(beta)*cos(beta),   -exp(1i*phi) + (exp(1i*phi)-1) * sin(beta)*sin(beta)];
G^(J+1)*psi 