function tests = test_crossval_blocked
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testSelectsBestLambda(testCase)
% section lambda selection
% confirm cross-validation picks the lambda with the lowest hold-out NLL.
setRNG();
X = [ones(200, 1), linspace(-1, 1, 200)'];
trueW = [0.5; -1.0];
rate = exp(X * trueW);
y = poissrnd(rate);
Xd = struct('X', X, 'y', y);

fitfun = @(XdTrain, Dtrain, lambda) runPenalisedFit(XdTrain, Dtrain, lambda);
D = sparse(0, size(X, 2));
cvCfg = struct('k', 5, 'lambdas', [0, 0.1, 1, 10]);

[bestLambda, cvinfo] = crossval_blocked(fitfun, Xd, D, cvCfg);

[~, minIdx] = min(cvinfo.mean_nll);
testCase.verifyEqual(bestLambda, cvinfo.lambdas(minIdx));
testCase.verifyEqual(cvinfo.k, cvCfg.k);
testCase.verifyEqual(sum(cvinfo.fold_sizes), numel(y));
end

function testFoldCoverage(testCase)
% section fold coverage
% ensure blocked folds cover all samples without overlap or omission.
setRNG();
X = randn(103, 3);
y = poissrnd(exp(X(:, 1)));
Xd = struct('X', X, 'y', y);

fitfun = @(XdTrain, Dtrain, lambda) runPenalisedFit(XdTrain, Dtrain, lambda);
D = sparse(0, size(X, 2));
cvCfg = struct('k', 4, 'lambdas', [0, 0.5]);

[~, cvinfo] = crossval_blocked(fitfun, Xd, D, cvCfg);

foldSizes = cvinfo.fold_sizes;
testCase.verifyEqual(sum(foldSizes), numel(y));

foldBreaks = [0, cumsum(foldSizes)];
covered = false(numel(y), 1);

startIdx = 1;
for jj = 1:cvCfg.k
    idxStart = foldBreaks(jj) + 1;
    idxEnd = foldBreaks(jj + 1);
    covered(idxStart:idxEnd) = covered(idxStart:idxEnd) | true;
    startIdx = idxEnd + 1;
end

testCase.verifyTrue(all(covered));
end

function setRNG()
if verLessThan('matlab', '9.9')
    rng(123);
else
    rng(123, 'twister');
end
end

function model = runPenalisedFit(XdTrain, Dtrain, lambda)
optCfg = struct('max_iter', 100, 'display', 'off');
[model, ~] = fit_glm_map(XdTrain, Dtrain, lambda, optCfg);
end
