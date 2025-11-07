function lfp_pred = predict_lfp(X, w, link)
%PREDICT_LFP Predict LFP signal from design matrix and weights.
%   lfp_pred = PREDICT_LFP(X, w, link) computes predicted LFP values
%   using the fitted GLM weights.
%
%   Inputs:
%     X    - design matrix (n_bins × n_features), sparse or dense
%     w    - weight vector (n_features × 1)
%     link - (optional) link function: 'identity' (default) or 'log'
%
%   Returns:
%     lfp_pred - predicted LFP signal (n_bins × 1)
%
%   Link functions:
%     'identity': lfp_pred = X * w (allows negative values)
%     'log':      lfp_pred = exp(X * w) (ensures positive values)

% section input validation and defaults
if nargin < 3 || isempty(link)
    link = 'identity';
end

w = w(:);

% section prediction based on link function
switch lower(link)
    case 'identity'
        % linear prediction: mu = X * w
        lfp_pred = X * w;

    case 'log'
        % exponential prediction: mu = exp(X * w)
        % clip linear predictor for numerical stability
        linpred = X * w;
        linpred = max(min(linpred, 50), -50);
        lfp_pred = exp(linpred);

    otherwise
        error('predict_lfp:InvalidLink', 'Link function must be "identity" or "log", got: %s', link);
end

end
