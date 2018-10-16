function betagal_refolding_kinetics(time, absorbance, a, rate, c, txtname)
    frac_absorbance = absorbance/100;    
    temp_x = time(1):1:time(end);
    params = [a rate c];
    [fitted_rate, r, J, covb, err, ErrorModel] = nlinfit(time, frac_absorbance, @kinetic_fit, params); 
    %ci = nlparci(fitparams, r, 'Jacobian', J); 
    betagal_abs_fit = kinetic_fit(fitted_rate, temp_x); 
    
    Fit.fit = betagal_abs_fit;
    Fit.params = fitted_rate;
    Fit.MSE = err;
    Fit.residuals = r;
    Fit.CovB = covb;
    StdDev = cov2corr(covb);
    Fit.StdDev = StdDev;
    Fit.Jacobian = J;
    Fit.ErrorModel = ErrorModel;
    
    figure
    plot(time, absorbance, 'ko');
    hold on
    plot(temp_x, betagal_abs_fit, 'k');
    
    eval(sprintf('%s = Fit', txtname));
    save(strcat(txtname, '.mat'),txtname);
end

function [ fit ] = kinetic_fit( params, time )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(time)    
    fit(i,1) = params(1)*(1 - exp(params(2)*time(i))) + params(3); 
end

end

