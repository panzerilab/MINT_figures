function [yFit, beta, R2, RMSE, asymp_val] = fit_info_curves(pop_size_ALL, mean_info, fit_function, info_scaling, SI_CI_info_ind)
%%
    % data to fit
    xData = pop_size_ALL;
    yData = mean_info;

    switch fit_function
        case 'quadratic'
            % Fit the data to a parabolic curve (quadratic function)
            p = polyfit(xData, yData, 2);
            % Generate a dense set of points from the fitted curve for plotting
            yFit = polyval(p, xData);  
            asymp_val = Nan;
            beta = 0;

        case 'logistic'
            % Logistic growth model function
            logisticModel = @(b,x) b(1) ./ (1 + exp(-b(2)*(x - b(3))));
            % Initial guess for parameters L, k, x0
            initialGuess = [max(yData), 1, median(xData)]; 
            % Fit model to data
            [beta,R,J,CovB,MSE] = nlinfit(xData, yData, logisticModel, initialGuess);
            asymp_val = beta(1);
            % Generate data for plotting the fitted curve
        %     xValues = linspace(min(xData), max(xData), 100);
            yFit = logisticModel(beta, xData);

        case 'Ince2013_model_2params'
            % Ince random model (fit two parameters)
            ince_model = @(phi, M) (1 - phi(1).^M) * phi(2);
            initialGuess = [0.1, 2];
            [beta, R, J, CovB, MSE] = nlinfit(xData, yData, ince_model, initialGuess);
            asymp_val = beta(2);
            yFit = ince_model(beta, xData);

        case 'Ince2013_model_1param'
            % Ince random model (fit one single parameter)
            if strcmp(info_scaling,'SI')
                phi2 = log2(2);
            elseif strcmp(info_scaling,'CI')
                phi2 = log2(2);
            else strcmp(info_scaling,'II')
                phi2 = nanmean(SI_CI_info_ind);
            end
            ince_model = @(phi, M) (1 - phi(1).^M) * phi2;
            initialGuess = [0.1];
            [beta, R, J, CovB, MSE] = nlinfit(xData, yData, ince_model, initialGuess);
            asymp_val = phi2;
            yFit = ince_model(beta, xData);            
    end

    % Goodness of fit calculations
    SStotal = sum((yData - mean(yData)).^2); % Total sum of squares
    SSresid = sum((yData - yFit).^2);    % Residual sum of squares
    R2 = 1 - SSresid/SStotal;       % Coefficient of determination
    RMSE = sqrt(mean((yData - yFit).^2)); % Root mean squared error
    
    % Display goodness of fit
    disp(['R-squared: ', num2str(R2)]);
    disp(['RMSE: ', num2str(RMSE)]);
    disp(['The asymptotic value (L) is: ', num2str(asymp_val)]);
end