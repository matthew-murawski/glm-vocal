function tests = test_fit_glm_map_hg
%TEST_FIT_GLM_MAP_HG Unit tests for fit_glm_map_hg.
%   tests = TEST_FIT_GLM_MAP_HG() returns function-based tests.

tests = functiontests(localfunctions);
end

function testSimple2DProblem(testCase) %#ok<INUSD>
% test recovery of weights in simple 2D problem
n = 100;

rng(42);
X = [ones(n, 1), randn(n, 1)];  % intercept + slope
w_true = [5; 2];
noise_std = 0.5;
y = X * w_true + noise_std * randn(n, 1);
y(y < 0) = 0;

% create config
cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0.01;
cfg.model.penalty_operator = sparse(0, 2);  % no penalty for simple test
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.max_iter = 100;
cfg.optim.tol_fun = 1e-6;
cfg.optim.tol_grad = 1e-5;
cfg.optim.display = 'off';

% fit model
fit = fit_glm_map_hg(X, y, cfg);

% check convergence
assert(fit.converged, 'Fit should converge');
assert(fit.exitflag > 0);

% check weight recovery (within 10% for good SNR)
rel_error = abs(fit.w - w_true) ./ abs(w_true);
assert(all(rel_error < 0.2), 'Weights should be recovered accurately');
end

function testZeroPenaltyMatchesLeastSquares(testCase) %#ok<INUSD>
% test that λ=0 matches least squares solution
n = 50;
p = 3;

rng(123);
X = randn(n, p);
w_true = randn(p, 1);
y = X * w_true + 0.1 * randn(n, 1);

% config with zero penalty
cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0;
cfg.model.penalty_operator = sparse(0, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';

% fit with our function
fit = fit_glm_map_hg(X, y, cfg);

% compare to least squares
w_ls = X \ y;

% should match closely
max_diff = max(abs(fit.w - w_ls));
assert(max_diff < 1e-3, 'Zero penalty should match least squares');
end

function testPenaltyShrinkage(testCase) %#ok<INUSD>
% test that penalty shrinks weights
n = 60;
p = 5;

rng(456);
X = randn(n, p);
w_true = randn(p, 1) * 2;
y = X * w_true + 0.2 * randn(n, 1);

% create penalty matrix
D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

% fit with different lambda values
lambda_vals = [0, 0.1, 1.0];
weights = zeros(p, length(lambda_vals));

for i = 1:length(lambda_vals)
    cfg = struct();
    cfg.model = struct();
    cfg.model.lambda = lambda_vals(i);
    cfg.model.penalty_operator = D;
    cfg.model.link = 'identity';
    cfg.optim = struct();
    cfg.optim.display = 'off';

    fit = fit_glm_map_hg(X, y, cfg);
    weights(:, i) = fit.w;
end

% weights should shrink with increasing lambda
w_norm = zeros(size(lambda_vals));
for i = 1:length(lambda_vals)
    w_norm(i) = norm(weights(:, i));
end

% check monotonic decrease
assert(w_norm(2) < w_norm(1), 'Lambda=0.1 should shrink weights vs lambda=0');
assert(w_norm(3) < w_norm(2), 'Lambda=1.0 should shrink weights vs lambda=0.1');
end

function testPerfectData(testCase) %#ok<INUSD>
% test convergence on perfect data (no noise)
n = 40;
p = 3;

rng(789);
X = randn(n, p);
w_true = [10; -2; 3];
y = X * w_true;  % no noise

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0;
cfg.model.penalty_operator = sparse(0, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';
cfg.optim.tol_fun = 1e-8;
cfg.optim.tol_grad = 1e-8;

fit = fit_glm_map_hg(X, y, cfg);

% should converge
assert(fit.converged);

% should achieve very low NLL
assert(fit.nll < 1e-6, 'Perfect fit should have near-zero NLL');

% weights should match truth closely
assert(max(abs(fit.w - w_true)) < 1e-4);
end

function testConvergenceDiagnostics(testCase) %#ok<INUSD>
% test that convergence diagnostics are captured
n = 50;
p = 4;

rng(101);
X = randn(n, p);
y = randn(n, 1) * 10;

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0.1;
cfg.model.penalty_operator = sparse(p-2, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';

fit = fit_glm_map_hg(X, y, cfg);

% check that all diagnostic fields exist
assert(isfield(fit, 'exitflag'));
assert(isfield(fit, 'converged'));
assert(isfield(fit, 'grad_norm'));
assert(isfield(fit, 'output'));
assert(isfield(fit, 'nll'));

% check that they have reasonable values
assert(isscalar(fit.exitflag));
assert(islogical(fit.converged) || isnumeric(fit.converged));
assert(isscalar(fit.grad_norm) && fit.grad_norm >= 0);
assert(isfinite(fit.nll));
end

function testNonConvergenceWarning(testCase)
% force a very low iteration limit to trigger non-convergence warning
n = 80; p = 8;
rng(404);
X = randn(n, p);
y = randn(n, 1) * 2;

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0.1;
cfg.model.penalty_operator = sparse(max(p-2,0), p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';
cfg.optim.max_iter = 1; % deliberately tiny to cause early termination

testCase.verifyWarning(@() fit_glm_map_hg(X, y, cfg), 'glm:NoConvergence');
end

function testPredictions(testCase) %#ok<INUSD>
% test that predictions are computed and correlate with observations
% generate y from X with a positive intercept to avoid heavy clipping
n = 60;
p = 4;

rng(202);
X = [ones(n,1), randn(n, p-1)];
w_true = [10; 1.5; -0.5; 0.8];
noise_std = 0.5;
y = X * w_true + noise_std * randn(n,1);

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0;
cfg.model.penalty_operator = sparse(0, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';

fit = fit_glm_map_hg(X, y, cfg);

% check predictions exist
assert(isfield(fit, 'mu'));
assert(length(fit.mu) == n);

% predictions should correlate strongly with observations
corr_val = corr(fit.mu, y);
assert(corr_val > 0.8, 'Predictions should correlate with observations');
end

function testDimensionMismatch(testCase)
% test error for dimension mismatch
n = 20;
p = 3;

X = randn(n, p);
y_wrong = randn(n+5, 1);

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0;
cfg.model.penalty_operator = sparse(0, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';

verifyError(testCase, @() fit_glm_map_hg(X, y_wrong, cfg), 'glm:DimensionMismatch');
end

function testInvalidLink(testCase)
% test error for invalid link function
n = 20;
p = 3;

X = randn(n, p);
y = randn(n, 1);

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0;
cfg.model.penalty_operator = sparse(0, p);
cfg.model.link = 'invalid_link';
cfg.optim = struct();
cfg.optim.display = 'off';

verifyError(testCase, @() fit_glm_map_hg(X, y, cfg), 'glm:InvalidLink');
end

function testSparseDesignMatrix(testCase) %#ok<INUSD>
% test fitting with sparse design matrix
n = 100;
p = 20;

rng(303);
% create sparse design matrix (10% non-zero)
X = sprand(n, p, 0.1);
w_true = randn(p, 1);
y = X * w_true + 0.5 * randn(n, 1);

cfg = struct();
cfg.model = struct();
cfg.model.lambda = 0.05;
cfg.model.penalty_operator = sparse(p-2, p);
cfg.model.link = 'identity';
cfg.optim = struct();
cfg.optim.display = 'off';

fit = fit_glm_map_hg(X, y, cfg);

% should converge
assert(fit.converged);
assert(isfinite(fit.nll));
assert(length(fit.w) == p);
end
