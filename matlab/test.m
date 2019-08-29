clc
close 
clear
A = [1/sqrt(2) 1/sqrt(2);
    1/sqrt(2) -1/sqrt(2)];
disp(A)
pretty(sym(A))
B = sqrt(A)
pretty(sym(B))