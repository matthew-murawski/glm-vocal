function tests = test_fit_glm_map
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testFitsSyntheticData(testCase)
% section synthetic fit
% verify the fitter recovers sensible parameters and positive rate predictions on synthetic data.
X = [ones(200, 1), linspace(-1, 1, 200)'];
trueW = [0.2; -0.4];
rate = exp(X * trueW);
setRNG();
y = poissrnd(rate);

Xd = struct('X', X, 'y', y);
D = sparse(0, size(X, 2));

lambda = 0.05;
optCfg = struct('max_iter', 200, 'tol_fun', 1e-8, 'tol_grad', 1e-7, 'display', 'off');

[wmap, fitinfo] = fit_glm_map(Xd, D, lambda, optCfg);

mu = exp(X * wmap.w);

initialNll = neglogli_poiss(zeros(size(wmap.w)), X, y, @(w) deal(0, zeros(size(w))));

testCase.verifyEqual(size(wmap.w), [size(X, 2), 1]);
testCase.verifyGreaterThan(all(mu), 0);
testCase.verifyGreaterThan(fitinfo.iterations, 0);
testCase.verifyLessThan(fitinfo.nll, initialNll);
end

function testPenalisedFitShrinksWeights(testCase)
% section penalty effect
% ensure adding a smoothness penalty shrinks weight magnitudes relative to unpenalised fit.
X = [ones(150, 1), linspace(-1, 1, 150)'];
wTrue = [0.3; 1.5];
rate = exp(X * wTrue);
setRNG();
y = poissrnd(rate);

Xd = struct('X', X, 'y', y);
D = sparse(1, 2, 1, 1, size(X, 2));

optCfg = struct('max_iter', 200, 'display', 'off');

[wFree, ~] = fit_glm_map(Xd, D, 0, optCfg);
[wPen, ~] = fit_glm_map(Xd, D, 5, optCfg);

penNorm = norm(wPen.w);
freeNorm = norm(wFree.w);

testCase.verifyLessThanOrEqual(penNorm, freeNorm);
end

function setRNG()
% section rng helper
% ensure deterministic sampling across test runs.
if verLessThan('matlab', '9.9')
    rng(42);
else
    rng(42, 'twister');
end
end
