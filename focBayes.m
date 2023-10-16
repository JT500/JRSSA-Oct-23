function [F,A,B] = focBayes(y,gamma,param,funH,funL)
% foc for 1-sided or 2-sided option problems
%
% INPUTS:
% y: value at which timing equation is evaluated
%
% beta: roots of fundamental quadratic. If beta is a scalar the function
% interprets beta<-.5 as a disinvestment options and beta>.5 as an investment
% option. If beta is a vector with one positive and one negative entry, a
% two-sided problem is assumed.
%
% funH: for one-sided problems this is the user-provided NPV and Jacobian.
% For two-sided problems it is the user-provided increasing payoff function
%
% funL: The user-provided decreasing payoff function in two-sided problems
% (not provided for one-sided problems)
%
% cycledum: dummy=0 if investment-abandonment decision, =1 if investment-mothballing decision 

xL = y(1);
xH = y(2);
g1 = gamma;
g2 = -gamma;

[Fnpv,Jnpv] = feval(funL,xL,param);
[FnpvH,JnpvH] = feval(funH,xH,param);

phiH = xL.^(.5+g1).*(1-xL).^(.5-g1);
JphiH = (xL./(1-xL)).^g1.*((.5+g1).*sqrt((1-xL)./xL)-(.5-g1).*sqrt(xL./(1-xL)));
phiL = xL.^(.5+g2).*(1-xL).^(.5-g2);
JphiL = (xL./(1-xL)).^g2.*((.5+g2).*sqrt((1-xL)./xL)-(.5-g2).*sqrt(xL./(1-xL)));

denom = phiL.*JphiH-JphiL.*phiH;
A = (phiL.*Jnpv-JphiL.*Fnpv)./denom;
B = (JphiH.*Fnpv-phiH.*Jnpv)./denom;

phiH1 = xH.^(.5+g1).*(1-xH).^(.5-g1);
JphiH1 = (xH./(1-xH)).^g1.*((.5+g1).*sqrt((1-xH)./xH)-(.5-g1).*sqrt(xH./(1-xH)));
phiL1 = xH.^(.5+g2).*(1-xH).^(.5-g2);
JphiL1 = (xH./(1-xH)).^g2.*((.5+g2).*sqrt((1-xH)./xH)-(.5-g2).*sqrt(xH./(1-xH)));

VL = A.* phiH1+B.*phiL1;
JVL = A.*JphiH1+B.*JphiL1;

F = [VL-FnpvH;JVL-JnpvH];


end
