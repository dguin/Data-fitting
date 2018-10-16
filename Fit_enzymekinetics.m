function Fit_enzymekinetics(Km, S0, Vmax, time, y, txtname)
% params(1) - Km, params(2) - S0, params(3) - Vmax

    fitparams0 = [Km, S0, Vmax];
    
    figure; set(gca, 'FontSize', 16); 
    plot(time, y, 'b', 'Linewidth', 2); 
    hold on; 
    [fitparams, r, J, covb, err, ErrorModel] = nlinfit(time, y, @EnzymeKinetics, fitparams0); 
    %ci = nlparci(fitparams, r, 'Jacobian', J); 

    Fit = EnzymeKinetics(fitparams, time); 

    Fit.fit = Fit;
    Fit.params = fitparams;
    Fit.MSE = err;
    Fit.residuals = r;
    Fit.CovB = covb;
    StdDev = cov2corr(covb);
    Fit.StdDev = StdDev;
    Fit.Jacobian = J;
    Fit.ErrorModel = ErrorModel;
    Fit.time = time;
        
    figure; set(gca, 'FontSize', 16); 
    plot(time, y, 'b', 'Linewidth', 2); 
    hold on; plot(Fit.time, Fit, 'k', 'Linewidth', 2);

    eval(sprintf('%s = Fit', txtname));
    save(strcat(txtname, '.mat'),txtname);
end



