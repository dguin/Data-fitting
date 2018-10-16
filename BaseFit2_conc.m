function BaseFit2_conc(varargin)

%thermx is the varied parameter, thermy is the thermodynamic
%data, lblimit is the fitting limit of the lower baseline, 
%hblimit is the fitting limit of the upper baseline

    switch nargin
        case 6 
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            Cm = cell2mat(varargin(3));
            lblimit = cell2mat(varargin(4));
            m1 = cell2mat(varargin(5));
            m2 = cell2mat(varargin(6));
            hblimit = 1;
            txtname = 'Fit';
            params_to_fit = [false false true true false false];
        case 7 
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            Cm = cell2mat(varargin(3));
            lblimit = cell2mat(varargin(4));
            m1 = cell2mat(varargin(5));
            m2 = cell2mat(varargin(6));
            hblimit = cell2mat(varargin(7));
            txtname = 'Fit';
            params_to_fit = [false false false false false false];
        case 8
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            Cm = cell2mat(varargin(3));
            lblimit = cell2mat(varargin(4));
            m1 = cell2mat(varargin(5));
            m2 = cell2mat(varargin(6));
            hblimit = cell2mat(varargin(7));
            txtname = cell2mat(varargin(8));
            params_to_fit = [false false false false false false];
        case 9
            thermx = cell2mat(varargin(1));
            thermy = cell2mat(varargin(2));
            Cm = cell2mat(varargin(3));
            lblimit = cell2mat(varargin(4));
            m1 = cell2mat(varargin(5));
            m2 = cell2mat(varargin(6));
            hblimit = cell2mat(varargin(7));
            txtname = cell2mat(varargin(8));
            params_to_fit = cell2mat(varargin(9));
        otherwise
            error('Error: not enough input arguements');
    end

    %estimate baselines 
    [lbase] = fit(thermx(1:lblimit), thermy(1:lblimit), 'poly1'); 
    [hbase] = fit(thermx(hblimit:end), thermy(hblimit:end), 'poly1'); 
    [m] = fit(thermx(m1:m2), thermy(m1:m2), 'poly1');
    %set up sigmoidal fit 

    concentration = [(thermx(1)):.025:(thermx(end))]'; 
    fitparams0 = [lbase.p1, lbase.p2, hbase.p1, hbase.p2, m.p1, Cm];

    figure; set(gca, 'FontSize', 16); 
    plot(thermx, thermy, 'bo', 'Linewidth', 2); 
    hold on; 
    [fitparams, r, J, covb, err, ErrorModel] = nlinfitsome(params_to_fit, thermx, thermy, @SigFit_conc, fitparams0); 
    %ci = nlparci(fitparams, r, 'Jacobian', J); 

    ConcFit = SigFit_conc(fitparams, concentration); 

    Fit.fit = ConcFit;
    Fit.params = fitparams;
    Fit.MSE = err;
    Fit.residuals = r;
    Fit.CovB = covb;
    StdDev = cov2corr(covb);
    Fit.StdDev = StdDev;
    Fit.Jacobian = J;
    Fit.ErrorModel = ErrorModel;
    Fit.conc = concentration;
    
    figure; set(gca, 'FontSize', 16); 
    plot(thermx, thermy, 'bo', 'Linewidth', 2); 
    hold on; plot(concentration, ConcFit, 'k', 'Linewidth', 2);

    eval(sprintf('%s = Fit', txtname));
    save(strcat(txtname, '.mat'),txtname);
end



