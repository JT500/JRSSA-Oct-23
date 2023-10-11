function [] = myplotbounds( piL, piH,mu,sig,prior,time )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% PLOT bounds and solution




lambdat=exp((mu/(sig^2))*(mu*1.1-(((mu)/2)*time)));
pit=((prior/(1-prior))*lambdat)/(1+(prior/(1-prior))*lambdat);

time_years=time;
posterior=pit;



X1=[0,1];
Y=[ piL piH];
YY=[ piL piH];

YMatrix=[Y; YY];

Y1=[piL piL];
Y2=[piH piH];

%plot(X1,Y1, 'k-');hold all;
%plot(X1,Y2, 'k-'); hold all;
%scatter(posterior,time_years);
%gridxy(0,piH);
%gridxy(0,piL);

createbounds50k(X1,YMatrix,time_years,posterior);
end

