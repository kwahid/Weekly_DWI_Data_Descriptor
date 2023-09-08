% Performs different types of fits on diffusion data for ADC
% fit_type:'ADC_exponential', 'ADC_linear_weighted', 'ADC_linear_simple', 'ADC_linear_fast'
function fit_output = fitParameter(parameter,fit_type,si,tr, userfile, ncoeffs, coeffs, tr_present,rsquared_threshold)

% input check
ok_ = isfinite(parameter) & isfinite(si);
if ~all( ok_ )
    warning( 'handle Nans And Infs by', ...
        'Ignoring NaNs and Infs in data' );
end

% check on data to avoid fitting of junk data
if(strcmp(fit_type,'ADC_linear_fast'))
    %Skip check as we are doing a fast linear fit
    r_squared = 2.0;

else
    ln_si = log(si);
    Ybar = mean(ln_si(ok_));
    Xbar = mean(parameter(ok_));
    y = ln_si(ok_)-Ybar;
    x = parameter(ok_)-Xbar;
    %     slope =sum(x.*y)/sum(x.^2);
    %     intercept = Ybar-slope.*Xbar; %#ok<NASGU>
    r_squared = (sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2)))^2;
    if ~isfinite(r_squared)
        r_squared = 0;
    end
end

% Continue on clean data
if r_squared>=rsquared_threshold
    if(strcmp(fit_type,'ADC_exponential'))
        % Restrict fits for ADC from 0 to Inf, and coefficient ('rho') from
        % 0 to inf
        fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 -Inf],'Upper',[Inf   0]);
        % The start point prevents convergance for some reason, do not use
        % 		st_ = [si(end) -.035 ];
        % 		set(fo_,'Startpoint',st_);
        %set(fo_,'Weight',w);
        ft_ = fittype('exp1');
        
        % Fit the model
        [cf_, gof] = fit(parameter(ok_),si(ok_),ft_,fo_);
        
        % Save Results
        sum_squared_error = gof.sse;
        r_squared = gof.rsquare;
        confidence_interval = confint(cf_,0.95);
        rho_fit = cf_.a;
        exponential_fit   = -1*cf_.b;
        exponential_95_ci = -1*confidence_interval(:,2);
    elseif(strcmp(fit_type,'ADC_linear_weighted'))
        % Restrict fits for ADC from 0 to Inf, and coefficient ('rho') from
        % 0 to inf
%         fo_ = fitoptions('method','LinearLeastSquares','Lower',[-1 -Inf],'Upper',[-Inf 0]);
        fo_ = fitoptions('method','LinearLeastSquares','Lower',[-1 -Inf],'Upper',[Inf 0]);
        ft_ = fittype('poly1');
        set(fo_,'Weight',si);
        ln_si = log(si);
        
        % Fit the model
        [cf_, gof] = fit(parameter(ok_),ln_si(ok_),ft_,fo_);
        
        % Save Results
        sum_squared_error = gof.sse;
        r_squared = gof.rsquare;
        confidence_interval = confint(cf_,0.95);
        rho_fit = cf_.p2;
        exponential_fit   = -1*cf_.p1;
        exponential_95_ci = -1*confidence_interval(:,1);
        
    elseif(strcmp(fit_type,'ADC_linear_simple'))
        % Restrict fits for T2 from 0 to Inf, and coefficient ('rho') from
        % 0 to inf
        fo_ = fitoptions('method','LinearLeastSquares','Lower',[-Inf 0],'Upper',[0 Inf]);
        ft_ = fittype('poly1');
        ln_si = log(si);
        
        % Fit the model
        [cf_, gof] = fit(parameter(ok_),ln_si(ok_),ft_,fo_);
        
        % Save Results
        sum_squared_error = gof.sse;
        r_squared = gof.rsquare;
        confidence_interval = confint(cf_,0.95);
        rho_fit = cf_.p2;
        exponential_fit   = -1*cf_.p1;
        exponential_95_ci = -1*confidence_interval(:,1);
    elseif(strcmp(fit_type,'ADC_linear_fast'))
        ln_si = log(si);
        
        % Fit the model
        Ybar = mean(ln_si(ok_));
        Xbar = mean(parameter(ok_));
        
        y = ln_si(ok_)-Ybar;
        x = parameter(ok_)-Xbar;
        slope =sum(x.*y)/sum(x.^2);
        intercept = Ybar-slope.*Xbar; 
        r_squared = (sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2)))^2;
        sum_squared_error = (1-r_squared)*sum(y.^2);
        if ~isfinite(r_squared) || ~isreal(r_squared)
            r_squared = 0;
        end
        % Save Results
        exponential_fit = -1*slope;
        rho_fit = intercept;
        % Confidence intervals not calculated
        exponential_95_ci(1) = -1;
        exponential_95_ci(2) = -1;
       
    elseif(strcmp(fit_type, 'user_input'))
        [PATHSTR,NAME,~] = fileparts(userfile);
        
        %Add the usefile function to path
%         path(path, PATHSTR)
%         
%         save('Moo.,mat', 'PATHSTR')
        
        userFN = str2func(NAME);
        
        % scale si, non-linear fit has trouble converging with big numbers
        scale_max = max(si);
        si = si./scale_max;
            
        [cf_, gof, output] = userFN(parameter(ok_), si(ok_));
        
        % Save Results
        sum_squared_error = gof.sse;
        r_squared = gof.rsquare;
        confidence_interval = confint(cf_,0.95);
        %rho_fit = cf_.a;
        %exponential_fit   = cf_.b;
        exponential_95_ci = confidence_interval;
      
        coeffvals = coeffvalues(cf_);
        coeffvals = coeffvals(:)';

    elseif(strcmp(fit_type,'none'))
        sum_squared_error = 0;
        r_squared = 1;
        rho_fit = 1;
        exponential_fit   = 1;
        exponential_95_ci = [1 1];
    end
    
else
    rho_fit = -2;
    exponential_fit = -2;
    exponential_95_ci(1) = -2;
    exponential_95_ci(2) = -2;
    sum_squared_error = -2;
    plus_c = -2;
    plus_c_low = -2;
    plus_c_high = -2;
end

if strcmp(fit_type, 'user_input')
    fit_output = [coeffvals r_squared exponential_95_ci(1,:) exponential_95_ci(2,:) sum_squared_error];
elseif strcmp(fit_type, 't2_exponential_plus_c')
    fit_output = [exponential_fit rho_fit r_squared exponential_95_ci(1) exponential_95_ci(2) sum_squared_error plus_c plus_c_low plus_c_high];
else
    fit_output = [exponential_fit rho_fit r_squared exponential_95_ci(1) exponential_95_ci(2) sum_squared_error];   
end
