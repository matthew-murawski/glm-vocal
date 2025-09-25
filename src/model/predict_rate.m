function mu = predict_rate(X, w)
% section poisson rate prediction
% compute the predicted firing rate mu = exp(X * w) with clipping for numerical stability.
w = w(:);
linpred = X * w;
linpred = max(min(linpred, 50), -50);
mu = exp(linpred);
end
