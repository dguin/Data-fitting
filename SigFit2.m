function [ThreeState] = SigFit2(params, tit)

%do a two state fit, where params are the fitting parameters
%and tit is the value varied (concentration or temperature)

%let deltaG1 = dG1*(T-Tm1), where deltaG1 = GL-GI 
%let deltaG2 = dG2*(T-Tm2), where deltaG2 = GI-GH
%fracLT = exp(-dG1*(T-Tm1)/kb*T)/(1+exp(-dG1*(T-Tm1)/kb*T)+exp(+dG2*(T-Tm2)/kb*T)) 
%fracHT = exp(+dG2*(T-Tm2)/kb*T)/(1+exp(-dG1*(T-Tm1)/kb*T)+exp(+dG2*(T-Tm2)/kb*T)) 
%fracIT = 1/(1+exp(-dG1*(T-Tm1)/kb*T)+exp(+dG2*(T-Tm2)/kb*T)) 


%LT baseline = lm*T+lb 
%HT baseline = hm*T+hb 
%IT baseline = im*T+ib

%S = fracLT*LTbase + fracHT*HTbase + fracIT*ITbase; 

%S = (lm*T+lb)*(1-fracHT-fracIT) + (hm*T+hb)*(fracHT) + (im*T+ib)*(fracIT)

%params = [lm, lb, hm, hb, im, ib, dG1, Tm1, dG2, Tm2] 

for i = 1:length(tit)
    LBase(i) = params(1)*tit(i) + params(2); 
    HBase1(i) = params(3)*tit(i) + params(4); 
    HBase2(i) = params(5)*tit(i) + params(6);
    K1(i) = exp(-params(7)*(tit(i)-params(8))/(8.31*tit(i)));
    K2(i) = exp(-params(9)*(tit(i)-params(10))/(8.31*tit(i)));
    
    ThreeState(i) = params(11)*(LBase(i)*(K1(i)/(1+K1(i))) + HBase1(i)*(1/(1+K1(i)))) + params(12)*(HBase1(i)*(K2(i)/(1+K2(i))) + HBase2(i)*(1/(1+K2(i)))); 

end

ThreeState = ThreeState';
plot(tit, ThreeState, 'k', 'Linewidth', 2); 
end
