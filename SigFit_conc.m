function [TwoState] = SigFit_conc(params, tit)

%do a two state fit, where params are the fitting parameters
%and tit is the value varied (concentration or temperature)

%let deltaG = dG1*(X-Xm), where deltaG = GH-GU 
%fracHX = exp(-dG1*(X-Xm)/kb*T)/(1+exp(-dG1*(X-Xm)/kb*T)) 

%LT baseline = lm*(X-Xm)+lb 
%HT baseline = hm*(X-Xm)+hb 

%S = fracLX*LXbase + fracHX*HXbase; 

%S = (lm*(X-Xm)+lb)*(1-fracHX) + (hm*(X-Xm)+hb)*(fracHX) 

%params = [lm, lb, hm, hb, dG1, Xm] 

for i = 1:length(tit)
    LBase(i) = params(1)*(tit(i)) + params(2); 
    HBase(i) = params(3)*(tit(i)) + params(4); 
    K(i) = exp((-params(5)*(params(6)-tit(i)))/(8.314*293));
    
    TwoState(i) = LBase(i)*(1/(1+K(i))) + HBase(i)*(K(i)/(1+K(i))); 

end

TwoState = TwoState';
plot(tit, TwoState, 'k', 'Linewidth', 2); 
end
