function BaseFit2(varargin)
%(thermx, thermy, dG1, Tm1, lblimit, hblimit, Tl, Th, tmelt, txtname, params_to_fit)
%thermx is the varied parameter, thermy is the thermodynamic
%data, lblimit is the fitting limit of the lower baseline, 
%hblimit is the fitting limit of the upper baseline
%for two state fits, set all unused parameters to 0

    switch nargin
        case 5 
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            dG1 = cell2mat(varargin(3));
            Tm1 = cell2mat(varargin(4));
            lblimit = cell2mat(varargin(5));
            hblimit = 1;
            tmelt = 1;
            txtname = 'Fit';
            params_to_fit = [false false true true false false true true];
        case 6 
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            dG1 = cell2mat(varargin(3));
            Tm1 = cell2mat(varargin(4));
            lblimit = cell2mat(varargin(5));
            hblimit = cell2mat(varargin(6));
            tmelt = 1;
            txtname = 'Fit';
            params_to_fit = [false false false false false false true true];
        case 8
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            dG1 = cell2mat(varargin(3));
            Tm1 = cell2mat(varargin(4));
            lblimit = cell2mat(varargin(5));
            hblimit = cell2mat(varargin(6));
            tmelt = cell2mat(varargin(7));
            txtname = cell2mat(varargin(8));
            params_to_fit = [false false false false false false true true];
        case 9
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            dG1 = cell2mat(varargin(3));
            Tm1 = cell2mat(varargin(4));
            lblimit = cell2mat(varargin(5));
            hblimit = cell2mat(varargin(6));
            tmelt = cell2mat(varargin(7));
            txtname = cell2mat(varargin(8));
            params_to_fit = cell2mat(varargin(9));
        otherwise
            error('Error: not enough input arguements');
    end

    %set up sigmoidal fit 
        
    Tl = thermx(1);
    Th = thermx(end);
    
    if tmelt ==1 
        thermx = thermx+273.15;
        Tm1 = Tm1+273.15;
        Tl = Tl+273.15;
        Th = Th+273.15;
    end 
    
    %estimate baselines 
    [lbase] = fit(thermx(1:lblimit)-(Tl), thermy(1:lblimit), 'poly1'); 
    [hbase] = fit(thermx(hblimit:end)-(Th), thermy(hblimit:end), 'poly1'); 

    tempfit = [round(thermx(1)):.5:round(thermx(end))]'; 

    %fitparams0 = [lbase.p1, lbase.p2, hbase.p1, 0, dG1, Tm1, Tl, Th];
    fitparams0 = [lbase.p1, lbase.p2, hbase.p1, hbase.p2, dG1, Tm1, Tl, Th];
    % fitparams0 = [lbase.p1, lbase.p2, hbase.p1, hbase.p2, dG1];

    %lb = [-Inf, -Inf, -Inf,  -Inf, -1000, thermx(1)]; 
    %hb = [Inf, Inf, Inf, Inf, 1000, thermx(end)]; 

    figure; set(gca, 'FontSize', 16); 
    plot(thermx, thermy, 'bo', 'Linewidth', 2); 
    hold on; 
    [fitparams, r, J, covb, err, ErrorModel] = nlinfitsome(params_to_fit, thermx, thermy, @SigFit, fitparams0); 
    %ci = nlparci(fitparams, r, 'Jacobian', J); 

    thermFit = SigFit(fitparams, tempfit); 

    Fit.fit = thermFit;
    Fit.params = fitparams;
    Fit.MSE = err;
    Fit.residuals = r;
    Fit.CovB = covb;
    StdDev = cov2corr(covb);
    Fit.StdDev = StdDev;
    Fit.Jacobian = J;
    Fit.ErrorModel = ErrorModel;
    Fit.params(6) = Fit.params(6) - 273.15;
    Fit.params(7:8) = Fit.params(7:8) - 273.15;
    %Fit.ci = ci;
    
    if tmelt == 1
        Fit.temps = tempfit-273.15;
        thermx = thermx-273.15;
    else
        Fit.temps = tempfit;
    end
        
    figure; set(gca, 'FontSize', 16); 
    plot(thermx, thermy, 'bo', 'Linewidth', 2); 
    hold on; plot(Fit.temps, thermFit, 'k', 'Linewidth', 2);

    eval(sprintf('%s = Fit', txtname));
    save(strcat(txtname, '.mat'),txtname);
end



