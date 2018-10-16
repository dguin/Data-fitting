function [Threestate] = SigFit3(params, tit)

%do a three state fit, where params are the fitting parameters
%and tit is the value varied (concentration or temperature)

%let deltaG1 = dG1*(T-Tm1), where deltaG1 = GL-GI -> K1
%let deltaG2 = dG2*(T-Tm2), where deltaG2 = GI-GH -> K2
%fracLT = 1/(1+K1+K1*K2)
%fracHT = K1*K2/(1+K1+K1*K2)
%fracIT = K1/(1+K1+K1*K2)


%LT baseline = lm*T+lb 
%HT baseline = hm*T+hb 
%IT baseline = (LT baseline + HT baseline)/2

%S = fracLT*LTbase + fracHT*HTbase + fracIT*ITbase; 
%See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2778058/

%params = [lm, lb, hm, hb, im, ib, dG1, Tm1, dG2, Tm2] 

for i = 1:length(tit)
    LBase(i) = params(1)*(tit(i)-params(11)) + params(2); 
    %HBase(i) = params(4);
    HBase(i) = params(3)*(tit(i)-params(12)) + params(4); 
    IBase(i) = params(5)*(tit(i)-params(12)-params(11)) + params(6); 
    K1(i) = (1/5)*exp(params(7)*(tit(i)-params(8))/(8.31*tit(i)));
    K2(i) = exp(params(9)*(tit(i)-params(10))/(8.31*tit(i)));
    %K1(i) = exp(-params(5)*(tit(i)-(273.15+35))/(8.31*tit(i)));
    
    Threestate(i) = LBase(i)*(1/(1+K1(i)+(K1(i)*K2(i)))) + HBase(i)*((K1(i)*K2(i))/(1+K1(i)+(K1(i)*K2(i)))) + IBase(i)*(K1(i)/(1+K1(i)+(K1(i)*K2(i)))) ; 
end

Threestate = Threestate';
plot(tit, Threestate, 'k', 'Linewidth', 2); 
end
