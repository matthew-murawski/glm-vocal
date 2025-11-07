function [nll, grad] = neglogli_gauss(w, X, y, L2op, link)
% NEGLOGLI_GAUSS Compute Gaussian negative log-likelihood and gradient.
%   [nll, grad] = NEGLOGLI_GAUSS(w, X, y, L2op, link) computes the negative
%   log-likelihood under a Gaussian (normal) distribution with optional
%   quadratic penalty and link function.
%
%   Inputs:
%     w     - weight vector (n_features × 1)
%     X     - design matrix (n_bins × n_features), sparse or dense
%     y     - response vector (n_bins × 1), continuous LFP values
%     L2op  - (optional) function handle for L2 penalty: [val, grad] = L2op(w)
%     link  - (optional) link function: 'identity' (default) or 'log'
%
%   Returns:
%     nll   - scalar negative log-likelihood (up to constant)
%     grad  - gradient vector (n_features × 1)
%
%   Model with identity link:
%     y ~ N(X*w, sigma^2)
%     NLL = (1/(2*sigma^2)) * ||y - X*w||^2
%     For fixed sigma, equivalent to minimizing sum of squared residuals
%
%   Model with log link:
%     y ~ N(exp(X*w), sigma^2)
%     Ensures positive predictions (useful for LFP power/envelope)

% section inputs and defaults
% compute the gaussian negative log-likelihood and gradient with optional quadratic penalty via L2op.
if nargin < 4 || isempty(L2op)
    L2op = @(w) deal(0, zeros(size(w)));
end
if nargin < 5 || isempty(link)
    link = 'identity';
end

if isempty(w)
    nll = 0;
    grad = zeros(0, 1);
    return
end

% ensure column vector inputs for compatibility
w = w(:);
y = y(:);

% compute linear predictor
if issparse(X)
    Xw = X * w;
else
    Xw = X * w;
end

% apply link function to get predicted mean
switch lower(link)
    case 'identity'
        % mu = X*w
        mu = Xw;
        % residuals: y - mu
        residuals = y - mu;
        % gradient contribution: -X' * residuals = X' * (mu - y)
        grad_data = -X' * residuals;

    case 'log'
        % mu = exp(X*w)
        % clip to avoid overflow
        Xw = max(min(Xw, 50), -50);
        mu = exp(Xw);
        % residuals: y - mu
        residuals = y - mu;
        % gradient: d/dw [0.5 * (y - exp(Xw))^2] = -(y - exp(Xw)) * exp(Xw) * X
        % = -residuals .* mu * X
        grad_data = -X' * (residuals .* mu);

    otherwise
        error('neglogli_gauss:InvalidLink', 'Link function must be "identity" or "log", got: %s', link);
end

% compute negative log-likelihood (up to constants)
% NLL = (1/(2*sigma^2)) * sum(residuals.^2)
% for optimization purposes, we can ignore 1/(2*sigma^2) scaling and constants
nll = 0.5 * sum(residuals.^2);

% gradient is already computed above
grad = grad_data;

% add penalty term
[penaltyVal, penaltyGrad] = L2op(w);
nll = nll + penaltyVal;
grad = grad + penaltyGrad;

end
