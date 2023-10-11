% XOLAIR case study 

clear all
clc

load('PSA_OUTPUT_INNOVATE')
QALY_SIM=deltaqaly;

% CLINICAL TRIAL COST ESTIMATE: 

% cost of running clinical trial 26000$ per patient - see
% http://www.prnewswire.com/news-releases/phase-3-clinical-trial-costs-exce
% ed-26000-per-patient-56447427.html

%ppc=(26000/1.57);  % cost per patient enrolled in $/Pound Sterling at 2012 exchange rate prices
%ppc=3600;
ppc=4200; % Computed by taking overall drug costs in Xolair and standard arms.( See Xolair modelvpopFeb2015)
r = .03;

patients = 9000; % population that would use xolair in England - source NICE report
g = 0;

cost0 = 26546;              % standard steroid treatment cost
cost1 = 67137;              % xolair cost
icost=cost1-cost0;          % Incremental cost

trial_pop=419;                  % total trial participants
%trial = 180/365;                  % trial length (26 weeks)
trial_cost = ppc*trial_pop;     % total trial cost

%ppc=(26000/1.57)*trial_pop/180;  % cost per patient enrolled in $/Pound Sterling at 2012 exchange rate prices
%ppc=(26000/1.57)/2;
%%% SET VALUE PER QALY
n=8;

for i=1:n

vQaly = 33500+i*2500;          % Value of QALY for NIMB>0 

MYQALY=QALY_SIM.*vQaly;
sig=std(MYQALY)*trial_pop*1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qaly1 = 11.37;          % Xolar QALY
Qaly0 = 10.17;          % standard care Qaly
iqaly=Qaly1-Qaly0;      % incremental qaly

%mu = (vQaly*(iqaly)-(icost))*trial_pop*2;     % mean NIMB
mu=((Qaly1-Qaly0)*vQaly-(cost1-cost0))*trial_pop*1e-6;
%mu=13223;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c =ppc*1e-6;                        % cost of adding a patient

% Payoffs
B1 =((Qaly1*vQaly-cost1)*patients/(r-g))*1e-6;

B0=((vQaly*Qaly0-cost1)*patients/(r-g))*1e-6;

B =((vQaly*Qaly0-cost0)*patients/(r-g))*1e-6;

ICER=(cost1-cost0)/(Qaly1-Qaly0);

I = 0;          %1*1e-6;
P = 0;
L = 0;

% Observation delay

delta = 0;
 
        g1 = .5*sqrt(1+4*r*(sig/mu)^2);
        g2 = -g1;

        %prior = .5;
        %zeta = prior/(1-prior);

        param = [B1;B0;B;P;L;I;c;r;delta];
        %param=param.*0.000001;
        if (B0-B>=I)||(B-L>=B1-I)||(B1-B0+P/r<=0)||(L-P/r<0)||(B+P/r<-(c/r-L))
            error('Inadmissible payoffs')
        end

        pibar = (B-B0+I)/(B1-B0+L);

        pi0 = [.1;.9];
        pi_ast = fsolve(@(y) focBayes(y,[g1;-g1],param,'HTA_Hprob','HTA_Lprob'),pi0);
        piL = pi_ast(1);
        piH = pi_ast(2);
       
        piL_lambda(i)=piL;
        piH_lambda(i)=piH;
        
        time=0.5;
        prior=0.5;
        lambdat=exp((mu/(sig^2))*(mu*1.1-(((mu)/2)*time)));
        pit=((prior/(1-prior))*lambdat)/(1+(prior/(1-prior))*lambdat);

        
        posterior(i)=pit;
end
% 
        % Display results
        disp('piL');
        disp(piL);
        disp('piH');
        disp(piH);
        
        
        
        
        
        
        x = 1:n;
        
        plot(piL_lambda)
        hold on
        plot(piH_lambda)
        hold on
        scatter(x,posterior)

%         
%         prior=0.5;
%         time=0.5;
% 
%         lambdat=exp((mu/(sig^2))*(mu*1.1-(((mu)/2)*time)));
%         pit=((prior/(1-prior))*lambdat)/(1+(prior/(1-prior))*lambdat);
%         
%         
%         B11 =((Qaly1*vQaly-cost1));
%         B00=((vQaly*Qaly0-cost1));
%               
%          
%         inv=piH*B11+(1-piH)*B00;
%         post_value=pit*B11+(1-pit)*B00;
%         inv_opt=inv-post_value;
%  
%         
%         ab=piL*(B00)+(1-piL)*B00;
%    
%         post_value_ab=pit*B11-(1-pit)*B00;
%         
%         
%         ab_option=post_value-ab;
%         
%         
%         % Payoff option
%         disp('Investment option')
%         disp(round(inv_opt))
%         
%         disp('Abandonment option')
%         disp(round(ab_option))
%         
%         disp('Investment payoff')
%         disp(round(inv))
%         
%         disp('Abandonment payoff')
%         disp(round(ab))
%         
%         disp('post_value_inv')
%         disp(round(post_value))
% 
%         
%         myplotbounds( piL, piH,mu,sig,prior,time )
%         
%         
%         %[V_pi] = value(piL,piH,prior,pit,B11,B00,P,r,sig,mu,c);
%         %myplotbounds( piL, piH,mu,sig )
%      