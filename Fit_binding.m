function [ output_args ] = Fit_binding (thermx, thermy, dG, m, txtname)
% %You’re titrating protein concentration P0 with variable dodine concentration 
% C. The 50:50 point in your curve shows up at dodine concentration C0. 
% So your measured data is a set of pairs (P0, C0).
% Now, we assume as a simplified model thateither 0 dodines or m dodines bind,
% and the protein is folded with 0, unfolded with m. So x is the amount of 
% bound (unfolded) protein. This is of course oversimplified since variable
% amounts of dodine may bind, and the protein may be partly unfolded!
% In the simple model, the folding-unfolding equilibrium is
% Keq= [(P0-x)(C-m*x)^m]/x  .   At the midpoint, C=C0 and x=P0/2, thus
% Keq = (C0-m*P0/2)^m = exp[-DG/RT]
% You have two adjustable parameters DG and m in that equation to make a 
% plot of P0 vs. C0 fit:
% P0= (2/m) * {C0-exp[-DG/(m*RT)]}, or making P0 the independent variable, 
% C0=P0*m/2+exp[-DG/(m*RT)]
% This is the uncooperative case, where binding of dodine neither enhances 
% nor diminishes further binding of dodine to the protein. When fitting it,
% you need to be very careful, as C0-exp[] can go negative if incorrect 
% values of DG and m are picked by the automatic fitting algorithm. 
% IGOR has ways of avoiding it, I assume Matlab can deal with it, 
% although I have no idea how, and Matlabs fitting algorithm is of inferior stability.
% If there is binding cooperativity, but m*P0/2 << C0 is still true, 
% we can assume a cooperative free energy, DG = DG0+dGc*C0. 
% Strictly speaking, C0 should be replaced by C0-m*P0/2 in that formula, 
% but then P0 is on both sides of the equation, and must be solved implicitly by iteration.) 
% Thus for (anti)cooperativity, the approximate fitting equation P0=f(C0) is
% P0= (2/m) * {C0-exp[-(DG0+dGc*C0)/(m*RT)]}
% with adjustable parameters m, DG0 and dGc. Again, care must be taken to 
% already pick initial conditions for the fit such that the resulting curve 
% is close to the experimental curve.

tempfit = [(thermx(1)):((thermx(end))-(thermx(1)))/50:(thermx(end))]';

fitparams0 = [dG, m];

[fitparams, r, J, covb, err] = nlinfit(thermx, thermy, @binding, fitparams0); 
ci = nlparci(fitparams, r, 'Jacobian', J); 

bindFit = binding(fitparams, tempfit); 

figure; set(gca, 'FontSize', 16); 
plot(thermx, thermy, 'bo', 'Linewidth', 2); 
hold on; plot(tempfit, bindFit, 'k', 'Linewidth', 2);

Fit.fit = bindFit;
Fit.params = fitparams;
Fit.ci = ci;
Fit.x = tempfit;

eval(sprintf('%s = Fit', txtname));
save(strcat(txtname, '_fit.mat'));


end

