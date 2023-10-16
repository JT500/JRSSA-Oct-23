function [F,J] = HTA_Lprob(y,param)

B1 = param(1);
B0 = param(2);
B = param(3);
P = param(4);
L = param(5);
I = param(6);
c = param(7);
r = param(8);
delta = param(9);

F = -y*(L-P/r)+c/r+B;
J = -(L-P/r)*ones(length(y),1);