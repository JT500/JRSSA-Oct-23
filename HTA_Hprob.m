function [F,J] = HTA_Hprob(y,param)

B1 = param(1);
B0 = param(2);
B = param(3);
P = param(4);
L = param(5);
I = param(6);
c = param(7);
r = param(8);
delta = param(9);

F = y*(B1-B0+P/r)+B0+c/r-I;
J = (B1-B0+P/r)*ones(length(y),1);