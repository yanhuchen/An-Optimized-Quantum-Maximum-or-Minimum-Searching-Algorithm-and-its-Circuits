clc
close
clear

q=15;
s=1;
a=0.999;
for i=0:q
    s= s*(1+a^(2^i));
end
s