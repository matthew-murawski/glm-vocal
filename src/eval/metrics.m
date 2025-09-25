function out = metrics(y, mu)
% section poisson metrics
% compute negative log-likelihood, pseudo-r2, and deviance explained relative to a constant-rate baseline.
y = y(:);
mu = mu(:);
if numel(y) ~= numel(mu)
    error('metrics:SizeMismatch', 'y and mu must have the same number of elements.');
end

mu = max(mu, eps);
baselineRate = max(mean(y), eps);

LL_const = gammaln(y + 1);
nll = sum(mu - y .* log(mu) + LL_const);

LL_model = sum(y .* log(mu) - mu);

mu_null = baselineRate * ones(size(y));
LL_null = sum(y .* log(mu_null) - mu_null);

positiveMask = y > 0;
LL_sat = sum(y(positiveMask) .* log(y(positiveMask)) - y(positiveMask));

denomPseudo = LL_sat - LL_null;
if denomPseudo == 0
    pseudoR2 = NaN;
else
    pseudoR2 = 1 - (LL_sat - LL_model) / denomPseudo;
end

dev_model = 2 * (LL_sat - LL_model);
dev_null = 2 * (LL_sat - LL_null);
if dev_null == 0
    devianceExplained = NaN;
else
    devianceExplained = 1 - (dev_model / dev_null);
end

out = struct();
out.nll = nll;
out.pseudoR2 = pseudoR2;
out.deviance_explained = devianceExplained;
end
