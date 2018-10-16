function [EnzKin] = EnzymeKinetics(params, time)

%Fits enzyme kinetics for the integrated function \frac{[S]}{K_\mathrm{M}} = W(F(t))\,  
%and F(t) = \frac{[S]_0}{K_\mathrm{M}}
%\exp\!\left(\frac{[S]_0}{K_\mathrm{M}} - \frac{V_\max}{K_\mathrm{M}}\,t
%\right) \,. From https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics#cite_note-27
%

EnzKin = zeros(length(time),1);

for i = 1:length(time)    
    phi = params(2)/params(1);
    EnzKin(i,1) = params(1)*lambertw(phi*exp(phi-(time(i)*(params(3)/params(1)))));
end

plot(time, EnzKin, 'k', 'Linewidth', 2); 
end
