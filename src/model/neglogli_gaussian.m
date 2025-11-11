function [nll, grad] = neglogli_gaussian(w, X, y, lambda, D)
%NEGLOGLI_GAUSSIAN Negative log-likelihood for Gaussian GLM (identity link).
%   [nll, grad] = NEGLOGLI_GAUSSIAN(w, X, y, lambda, D) computes the
%   negative log-likelihood and gradient for a Gaussian GLM with identity
%   link function and L2 smoothness penalty.
%
%   Mathematical model:
%       Linear predictor: η = X * w
%       Mean: μ = η (identity link)
%       Likelihood: y ~ Normal(μ, σ²)
%       NLL = (1/2n) * sum((y - μ)²) + (λ/2) * sum((D*w)²)
%
%   Inputs:
%       w       - [p × 1] weight vector
%       X       - [n × p] design matrix (can be sparse)
%       y       - [n × 1] observed power trace
%       lambda  - scalar regularization strength (≥ 0)
%       D       - penalty operator (second-difference matrix)
%
%   Outputs:
%       nll     - scalar negative log-likelihood
%       grad    - [p × 1] gradient vector
%
%   Computation:
%       Residuals: r = y - X*w
%       NLL: (1/(2*n)) * sum(r²) + (λ/2) * sum((D*w)²)
%       Gradient: -(1/n) * X' * r + λ * D' * D * w

% validate inputs
validate_glm_inputs(w, X, y, lambda, D);

% ensure column vectors
w = w(:);
y = y(:);

n = length(y);

% compute linear predictor and residuals
% use sparse matrix multiplication if X is sparse
mu = X * w;
r = y - mu;

% compute negative log-likelihood (data term)
% normalize by n for scale invariance
nll_data = sum(r.^2) / (2 * n);

% compute penalty term
% use matrix-vector product D*w rather than forming D'*D explicitly
Dw = D * w;
nll_penalty = (lambda / 2) * sum(Dw.^2);

% total negative log-likelihood
nll = nll_data + nll_penalty;

% compute gradient if requested
if nargout > 1
    % data term: -(1/n) * X' * r
    % this is the negative gradient of sum((y - X*w)²) / (2n)
    grad_data = -(X' * r) / n;

    % penalty term: λ * D' * D * w
    % use matrix-vector products to avoid forming D'*D
    grad_penalty = lambda * (D' * Dw);

    % total gradient
    grad = grad_data + grad_penalty;
end

end

%% helper functions

function validate_glm_inputs(w, X, y, lambda, D)
%VALIDATE_GLM_INPUTS Check dimensions and validity of GLM inputs.

% check that inputs are provided
if isempty(w)
    error('glm:InvalidInput', 'Weight vector w cannot be empty.');
end
if isempty(X)
    error('glm:InvalidInput', 'Design matrix X cannot be empty.');
end
if isempty(y)
    error('glm:InvalidInput', 'Response vector y cannot be empty.');
end

% ensure w is a vector
if ~isvector(w)
    error('glm:InvalidInput', 'Weight vector w must be a vector, got size [%s].', ...
        num2str(size(w)));
end

% ensure y is a vector
if ~isvector(y)
    error('glm:InvalidInput', 'Response vector y must be a vector, got size [%s].', ...
        num2str(size(y)));
end

% get dimensions
p = length(w);
n = length(y);

% check X dimensions
[n_rows, n_cols] = size(X);
if n_rows ~= n
    error('glm:DimensionMismatch', ...
        'Design matrix X has %d rows but y has %d elements.', n_rows, n);
end
if n_cols ~= p
    error('glm:DimensionMismatch', ...
        'Design matrix X has %d columns but w has %d elements.', n_cols, p);
end

% check that y contains only finite values
if any(~isfinite(y))
    error('glm:InvalidInput', 'Response vector y contains non-finite values.');
end

% check lambda
if ~isscalar(lambda)
    error('glm:InvalidInput', 'Regularization strength lambda must be a scalar.');
end
if ~isfinite(lambda)
    error('glm:InvalidInput', 'Regularization strength lambda must be finite.');
end
if lambda < 0
    error('glm:InvalidInput', 'Regularization strength lambda must be non-negative, got %.4f.', lambda);
end

% check D dimensions if penalty is being applied
if lambda > 0 && ~isempty(D)
    [d_rows, d_cols] = size(D);
    if d_cols ~= p
        error('glm:DimensionMismatch', ...
            'Penalty operator D has %d columns but w has %d elements.', d_cols, p);
    end
end

end
