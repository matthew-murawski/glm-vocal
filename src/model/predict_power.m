function mu = predict_power(w, X, link)
%PREDICT_POWER Predict mean power from GLM weights.
%   mu = PREDICT_POWER(w, X, link) computes predicted mean power using
%   the specified link function.
%
%   Inputs:
%       w    - [p × 1] weight vector
%       X    - [n × p] design matrix
%       link - 'identity' or 'log'
%
%   Output:
%       mu   - [n × 1] predicted mean power (always non-negative)
%
%   Link functions:
%       'identity': mu = X * w
%       'log':      mu = exp(X * w)

% validate inputs
if ~isvector(w)
    error('glm:InvalidInput', 'w must be a vector.');
end
if ~ismatrix(X)
    error('glm:InvalidInput', 'X must be a matrix.');
end

% ensure w is column vector
w = w(:);

% get dimensions
[n, p] = size(X);
if length(w) ~= p
    error('glm:DimensionMismatch', ...
        'X has %d columns but w has %d elements.', p, length(w));
end

% compute linear predictor
eta = X * w;

% apply link function
switch lower(link)
    case 'identity'
        mu = eta;

    case 'log'
        % clip eta to prevent overflow in exp
        eta_clipped = min(eta, 20);  % exp(20) ≈ 5e8
        if any(eta > 20)
            warning('glm:LinkOverflow', ...
                'Linear predictor exceeds 20 in %d bins. Clipping to prevent overflow.', ...
                sum(eta > 20));
        end
        mu = exp(eta_clipped);

    otherwise
        error('glm:InvalidLink', ...
            'Link function must be ''identity'' or ''log''. Got: %s', link);
end

% ensure non-negative predictions
if any(mu < 0)
    n_negative = sum(mu < 0);
    warning('glm:NegativePredictions', ...
        'Predicted power is negative in %d bins (%.2f%%). Setting to zero.', ...
        n_negative, 100 * n_negative / n);
    mu(mu < 0) = 0;
end

end
