function [partial_resid, fitted] = get_partial_residuals(model, data, predictor)
    % Get model formula
    formula = model.Formula;
    
    % Fit model without predictor of interest
    formula_without = regexprep(formula, [predictor ' \+ '], '');
    model_without = fitlme(data, formula_without);
    
    % Calculate partial residuals
    partial_resid = model.response - (model_without.fitted - model_without.offset);
    fitted = data.(predictor) * model.Coefficients.Estimate(strcmp(model.Coefficients.Name, predictor));
end