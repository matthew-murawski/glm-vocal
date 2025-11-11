function tests = test_neglogli_gaussian
%TEST_NEGLOGLI_GAUSSIAN Unit tests for neglogli_gaussian.
%   tests = TEST_NEGLOGLI_GAUSSIAN() returns function-based tests.

tests = functiontests(localfunctions);
end

function testNoPenalty(testCase) %#ok<INUSD>
% test that λ=0 matches manual sum of squared residuals
n = 50;
p = 3;

rng(42);
X = randn(n, p);
w = randn(p, 1);
y = X * w + 0.1 * randn(n, 1);

lambda = 0;
D = sparse(p-2, p);  % dummy penalty matrix

[nll, grad] = neglogli_gaussian(w, X, y, lambda, D);

% manual computation
r = y - X * w;
nll_expected = sum(r.^2) / (2 * n);

assert(abs(nll - nll_expected) < 1e-10);
end

function testGradientNoPenalty(testCase) %#ok<INUSD>
% test gradient against numerical gradient when λ=0
n = 30;
p = 4;

rng(123);
X = randn(n, p);
w = randn(p, 1);
y = X * w + 0.2 * randn(n, 1);

lambda = 0;
D = sparse(p-2, p);

[~, grad_analytical] = neglogli_gaussian(w, X, y, lambda, D);

% numerical gradient via finite differences
eps = 1e-7;
grad_numerical = zeros(p, 1);
for i = 1:p
    w_plus = w;
    w_plus(i) = w_plus(i) + eps;
    nll_plus = neglogli_gaussian(w_plus, X, y, lambda, D);

    w_minus = w;
    w_minus(i) = w_minus(i) - eps;
    nll_minus = neglogli_gaussian(w_minus, X, y, lambda, D);

    grad_numerical(i) = (nll_plus - nll_minus) / (2 * eps);
end

% check relative error
rel_error = max(abs(grad_analytical - grad_numerical) ./ (abs(grad_numerical) + 1e-10));
assert(rel_error < 1e-4);
end

function testWithPenalty(testCase) %#ok<INUSD>
% test that penalty increases NLL
n = 40;
p = 5;

rng(456);
X = randn(n, p);
w = randn(p, 1);
y = X * w + 0.15 * randn(n, 1);

% create second-difference matrix
D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

% compute NLL with and without penalty
lambda0 = 0;
nll0 = neglogli_gaussian(w, X, y, lambda0, D);

lambda1 = 1.0;
nll1 = neglogli_gaussian(w, X, y, lambda1, D);

% NLL should increase with penalty
assert(nll1 > nll0);
end

function testGradientWithPenalty(testCase) %#ok<INUSD>
% test gradient against numerical gradient with penalty
n = 25;
p = 6;

rng(789);
X = randn(n, p);
w = randn(p, 1);
y = X * w + 0.1 * randn(n, 1);

% create second-difference matrix
D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

lambda = 0.5;

[~, grad_analytical] = neglogli_gaussian(w, X, y, lambda, D);

% numerical gradient
eps = 1e-7;
grad_numerical = zeros(p, 1);
for i = 1:p
    w_plus = w;
    w_plus(i) = w_plus(i) + eps;
    nll_plus = neglogli_gaussian(w_plus, X, y, lambda, D);

    w_minus = w;
    w_minus(i) = w_minus(i) - eps;
    nll_minus = neglogli_gaussian(w_minus, X, y, lambda, D);

    grad_numerical(i) = (nll_plus - nll_minus) / (2 * eps);
end

% check relative error
rel_error = max(abs(grad_analytical - grad_numerical) ./ (abs(grad_numerical) + 1e-10));
assert(rel_error < 1e-4);
end

function testGradientAtMultiplePoints(testCase) %#ok<INUSD>
% test gradient accuracy at multiple random w values
n = 30;
p = 4;

rng(101);
X = randn(n, p);
y = randn(n, 1);

D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

lambda = 0.1;
eps = 1e-7;

% test at 5 random points
for trial = 1:5
    w = randn(p, 1);
    [~, grad_analytical] = neglogli_gaussian(w, X, y, lambda, D);

    % numerical gradient
    grad_numerical = zeros(p, 1);
    for i = 1:p
        w_plus = w;
        w_plus(i) = w_plus(i) + eps;
        nll_plus = neglogli_gaussian(w_plus, X, y, lambda, D);

        w_minus = w;
        w_minus(i) = w_minus(i) - eps;
        nll_minus = neglogli_gaussian(w_minus, X, y, lambda, D);

        grad_numerical(i) = (nll_plus - nll_minus) / (2 * eps);
    end

    rel_error = max(abs(grad_analytical - grad_numerical) ./ (abs(grad_numerical) + 1e-10));
    assert(rel_error < 1e-4);
end
end

function testAllZeros(testCase) %#ok<INUSD>
% test edge case: all zeros
n = 20;
p = 3;

X = zeros(n, p);
y = zeros(n, 1);
w = zeros(p, 1);

lambda = 0;
D = sparse(p-2, p);

[nll, grad] = neglogli_gaussian(w, X, y, lambda, D);

% should be zero
assert(nll == 0);
assert(all(grad == 0));
end

function testPerfectFit(testCase) %#ok<INUSD>
% test edge case: perfect fit with no noise
n = 50;
p = 3;

rng(222);
X = randn(n, p);
w_true = randn(p, 1);
y = X * w_true;  % no noise

lambda = 0;
D = sparse(p-2, p);

[nll, grad] = neglogli_gaussian(w_true, X, y, lambda, D);

% NLL should be zero (perfect fit)
assert(abs(nll) < 1e-10);
% gradient should be zero at minimum
assert(max(abs(grad)) < 1e-10);
end

function testDimensionMismatchError(testCase)
% test that dimension mismatches trigger errors
n = 20;
p = 5;

X = randn(n, p);
y = randn(n, 1);
w_wrong = randn(p+2, 1);  % wrong size

lambda = 0;
D = sparse(p-2, p);

verifyError(testCase, @() neglogli_gaussian(w_wrong, X, y, lambda, D), 'glm:DimensionMismatch');
end

function testYDimensionMismatchError(testCase)
% test that y dimension mismatch triggers error
n = 20;
p = 5;

X = randn(n, p);
y_wrong = randn(n+5, 1);  % wrong size
w = randn(p, 1);

lambda = 0;
D = sparse(p-2, p);

verifyError(testCase, @() neglogli_gaussian(w, X, y_wrong, lambda, D), 'glm:DimensionMismatch');
end

function testSparseVsDense(testCase) %#ok<INUSD>
% test that sparse and dense X produce same results
n = 40;
p = 8;

rng(333);
X_dense = randn(n, p);
X_sparse = sparse(X_dense);

w = randn(p, 1);
y = X_dense * w + 0.1 * randn(n, 1);

D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

lambda = 0.2;

[nll_dense, grad_dense] = neglogli_gaussian(w, X_dense, y, lambda, D);
[nll_sparse, grad_sparse] = neglogli_gaussian(w, X_sparse, y, lambda, D);

assert(abs(nll_dense - nll_sparse) < 1e-10);
assert(max(abs(grad_dense - grad_sparse)) < 1e-10);
end

function testLargeSparseMatrix(testCase) %#ok<INUSD>
% test memory efficiency with large sparse matrix
n = 1000;
p = 50;

% create sparse design matrix (only 5% non-zero)
rng(444);
X = sprand(n, p, 0.05);

w = randn(p, 1);
y = X * w + 0.1 * randn(n, 1);

D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

lambda = 0.1;

% should run without memory issues
[nll, grad] = neglogli_gaussian(w, X, y, lambda, D);

assert(isfinite(nll));
assert(all(isfinite(grad)));
assert(length(grad) == p);
end

function testNegativeLambdaError(testCase)
% test that negative lambda triggers error
n = 20;
p = 3;

X = randn(n, p);
y = randn(n, 1);
w = randn(p, 1);

lambda = -0.5;  % invalid
D = sparse(p-2, p);

verifyError(testCase, @() neglogli_gaussian(w, X, y, lambda, D), 'glm:InvalidInput');
end

function testNonFiniteYError(testCase)
% test that non-finite y triggers error
n = 20;
p = 3;

X = randn(n, p);
y = randn(n, 1);
y(5) = NaN;  % invalid
w = randn(p, 1);

lambda = 0;
D = sparse(p-2, p);

verifyError(testCase, @() neglogli_gaussian(w, X, y, lambda, D), 'glm:InvalidInput');
end

function testEmptyInputsError(testCase)
% test that empty inputs trigger errors
n = 20;
p = 3;

X = randn(n, p);
y = randn(n, 1);

lambda = 0;
D = sparse(p-2, p);

% empty w
verifyError(testCase, @() neglogli_gaussian([], X, y, lambda, D), 'glm:InvalidInput');

% empty X
verifyError(testCase, @() neglogli_gaussian(randn(p,1), [], y, lambda, D), 'glm:InvalidInput');

% empty y
verifyError(testCase, @() neglogli_gaussian(randn(p,1), X, [], lambda, D), 'glm:InvalidInput');
end
