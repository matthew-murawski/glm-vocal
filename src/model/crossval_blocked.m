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

n = numel(y);
foldSizes = computeFoldSizes(n, k);

foldBreaks = [0, cumsum(foldSizes)];

nLambda = numel(lambdas);
freqNll = zeros(nLambda, k);

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

        if isempty(D)
            Dtrain = D;
        else
            Dtrain = D;
        end

        model = fitfun(XdTrain, Dtrain, lambda);
        w = model.w;

        Xw = Xtest * w;
        Xw = max(min(Xw, 50), -50);
        mu = exp(Xw);
        nll = sum(mu - ytest .* Xw);
        freqNll(ll, fold) = nll / numel(ytest);
    end
end

meanNll = mean(freqNll, 2);
[~, bestIdx] = min(meanNll);
best_lambda = lambdas(bestIdx);

cvinfo = struct();
cvinfo.lambdas = lambdas;
cvinfo.mean_nll = meanNll;
cvinfo.fold_nll = freqNll;
cvinfo.k = k;
cvinfo.fold_sizes = foldSizes;
end

function foldSizes = computeFoldSizes(n, k)
% section fold helper
% compute sizes of contiguous folds covering all samples.
baseSize = floor(n / k);
remainder = mod(n, k);
foldSizes = repmat(baseSize, 1, k);
foldSizes(1:remainder) = foldSizes(1:remainder) + 1;
end
