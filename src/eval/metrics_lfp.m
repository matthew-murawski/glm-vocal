function out = metrics_lfp(y, lfp_pred)
%METRICS_LFP Compute evaluation metrics for continuous LFP predictions.
%   out = METRICS_LFP(y, lfp_pred) computes mean squared error, R-squared,
%   correlation, and other metrics for continuous signal prediction.
%
%   Inputs:
%     y        - actual LFP signal (n_bins × 1)
%     lfp_pred - predicted LFP signal (n_bins × 1)
%
%   Returns:
%     out - struct with fields:
%           .mse               - mean squared error
%           .rmse              - root mean squared error
%           .r2                - R-squared (fraction of variance explained)
%           .correlation       - Pearson correlation coefficient
%           .mae               - mean absolute error
%           .explained_variance - explained variance score

% section input validation
y = y(:);
lfp_pred = lfp_pred(:);
if numel(y) ~= numel(lfp_pred)
    error('metrics_lfp:SizeMismatch', 'y and lfp_pred must have the same number of elements.');
end

n = numel(y);
if n == 0
    % handle empty case
    out = struct('mse', NaN, 'rmse', NaN, 'r2', NaN, 'correlation', NaN, ...
                 'mae', NaN, 'explained_variance', NaN);
    return
end

% section residuals and error metrics
residuals = y - lfp_pred;

% mean squared error
mse = mean(residuals.^2);

% root mean squared error
rmse = sqrt(mse);

% mean absolute error
mae = mean(abs(residuals));

% section variance-based metrics
% R-squared: fraction of variance explained
y_mean = mean(y);
SS_tot = sum((y - y_mean).^2);
SS_res = sum(residuals.^2);

if SS_tot == 0
    % constant signal: R² is undefined
    r2 = NaN;
else
    r2 = 1 - (SS_res / SS_tot);
end

% explained variance score (similar to R² but uses variance instead of sum of squares)
var_y = var(y);
var_residuals = var(residuals);

if var_y == 0
    explained_variance = NaN;
else
    explained_variance = 1 - (var_residuals / var_y);
end

% section correlation
% Pearson correlation coefficient between actual and predicted
if n > 1
    corr_matrix = corrcoef(y, lfp_pred);
    correlation = corr_matrix(1, 2);
else
    correlation = NaN;
end

% handle edge case where correlation is NaN due to constant signal
if isnan(correlation) && (var(y) == 0 || var(lfp_pred) == 0)
    % one or both signals are constant
    if var(y) == 0 && var(lfp_pred) == 0 && y_mean == mean(lfp_pred)
        correlation = 1;  % perfect agreement
    else
        correlation = 0;  % no correlation
    end
end

% section outputs
out = struct();
out.mse = mse;
out.rmse = rmse;
out.r2 = r2;
out.correlation = correlation;
out.mae = mae;
out.explained_variance = explained_variance;

end
