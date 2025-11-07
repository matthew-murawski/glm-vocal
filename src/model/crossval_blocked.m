function [best_lambda, cvinfo] = crossval_blocked(fitfun, Xd, D, cvCfg)
% section inputs and defaults
% perform blocked cross-validation over candidate lambdas, returning the best lambda and per-fold metrics.
if nargin < 4 || isempty(cvCfg)
    error('crossval_blocked:MissingConfig', 'cvCfg with fields k and lambdas is required.');
end

if ~isfield(cvCfg, 'k') || ~isfield(cvCfg, 'lambdas')
    error('crossval_blocked:MissingFields', 'cvCfg must include k and lambdas.');
end

k = double(cvCfg.k);
if ~isscalar(k) || k < 2
    error('crossval_blocked:InvalidK', 'k must be an integer >= 2.');
end

lambdas = double(cvCfg.lambdas(:)');
if isempty(lambdas)
    error('crossval_blocked:NoLambdas', 'lambdas array must be non-empty.');
end

X = Xd.X;
y = Xd.y;

% detect response type from Xd (default to 'spikes' for backward compatibility)
if isfield(Xd, 'response_type')
    response_type = Xd.response_type;
else
    response_type = 'spikes';
end

% determine evaluation metric based on response type
if strcmpi(response_type, 'lfp')
    metric_name = 'mse';  % mean squared error for continuous LFP
else
    metric_name = 'nll';  % negative log-likelihood for spike counts
end

n = numel(y);
foldSizes = computeFoldSizes(n, k);

foldBreaks = [0, cumsum(foldSizes)];

nLambda = numel(lambdas);
freqMetric = zeros(nLambda, k);

for ll = 1:nLambda
    lambda = lambdas(ll);
    for fold = 1:k
        idxStart = foldBreaks(fold) + 1;
        idxEnd = foldBreaks(fold + 1);
        testIdx = idxStart:idxEnd;
        trainIdx = setdiff(1:n, testIdx);

        Xtrain = X(trainIdx, :);
        ytrain = y(trainIdx);
        Xtest = X(testIdx, :);
        ytest = y(testIdx);

        XdTrain = struct('X', Xtrain, 'y', ytrain);
        % preserve response_type in training data struct
        if isfield(Xd, 'response_type')
            XdTrain.response_type = Xd.response_type;
        end

        if isempty(D)
            Dtrain = D;
        else
            Dtrain = D;
        end

        model = fitfun(XdTrain, Dtrain, lambda);
        w = model.w;

        % compute evaluation metric based on response type
        if strcmpi(response_type, 'lfp')
            % for LFP: use mean squared error
            % assume identity link (can extend to support log link if needed)
            lfp_pred = Xtest * w;
            residuals = ytest - lfp_pred;
            mse = sum(residuals.^2) / numel(ytest);
            freqMetric(ll, fold) = mse;
        else
            % for spikes: use Poisson negative log-likelihood
            Xw = Xtest * w;
            Xw = max(min(Xw, 50), -50);
            mu = exp(Xw);
            nll = sum(mu - ytest .* Xw);
            freqMetric(ll, fold) = nll / numel(ytest);
        end
    end
end

meanMetric = mean(freqMetric, 2);
[~, bestIdx] = min(meanMetric);
best_lambda = lambdas(bestIdx);

cvinfo = struct();
cvinfo.lambdas = lambdas;
cvinfo.metric_name = metric_name;
cvinfo.mean_metric = meanMetric;
cvinfo.fold_metric = freqMetric;
cvinfo.k = k;
cvinfo.fold_sizes = foldSizes;
cvinfo.response_type = response_type;

% for backward compatibility, also store with old field names
if strcmpi(response_type, 'spikes')
    cvinfo.mean_nll = meanMetric;
    cvinfo.fold_nll = freqMetric;
end
end

function foldSizes = computeFoldSizes(n, k)
% section fold helper
% compute sizes of contiguous folds covering all samples.
baseSize = floor(n / k);
remainder = mod(n, k);
foldSizes = repmat(baseSize, 1, k);
foldSizes(1:remainder) = foldSizes(1:remainder) + 1;
end
