function [nll, grad] = neglogli_poiss(w, X, y, L2op)
% section inputs and defaults
% compute the poisson negative log-likelihood and gradient with optional quadratic penalty via L2op.
if nargin < 4 || isempty(L2op)
    L2op = @(w) deal(0, zeros(size(w)));
end

if isempty(w)
    nll = 0;
    grad = zeros(0, 1);
    return
end

% ensure column vector inputs for compatibility
w = w(:);
y = y(:);

if issparse(X)
    Xw = X * w;
else
    Xw = X * w;
end

% clip the linear predictor to avoid overflow in the exponential
Xw = max(min(Xw, 50), -50);
mu = exp(Xw);

% compute the negative log-likelihood (up to additive constant independent of w)
nll = sum(mu - y .* Xw);
grad = X' * (mu - y);

[penaltyVal, penaltyGrad] = L2op(w);
nll = nll + penaltyVal;
grad = grad + penaltyGrad;
end
