function tests = test_neglogli_poiss
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testGradientMatchesFiniteDifference(testCase)
% section gradient check
% verify analytic gradient matches finite-difference approximation on a small synthetic problem.
X = [1 0; 1 1; 1 2];
y = [0; 1; 2];
w = [0.2; -0.1];

lambda = 0.5;
L2op = @(w) deal(lambda * (w' * w), 2 * lambda * w);

[nll, grad] = neglogli_poiss(w, X, y, L2op);
	epsVal = 1e-6;
finiteGrad = zeros(size(w));
for ii = 1:numel(w)
    dw = zeros(size(w));
    dw(ii) = epsVal;
    wPlus = w + dw;
    wMinus = w - dw;
    nllPlus = neglogli_poiss(wPlus, X, y, L2op);
    nllMinus = neglogli_poiss(wMinus, X, y, L2op);
    finiteGrad(ii) = (nllPlus - nllMinus) / (2 * epsVal);
end

testCase.verifyEqual(nll, neglogli_poiss(w, X, y, L2op));
testCase.verifyLessThan(norm(grad - finiteGrad), 1e-4);
end

function testHandlesSparseInput(testCase)
% section sparse handling
% ensure the implementation works with sparse design matrices.
X = sparse([1 0 0; 1 1 0; 0 0 1]);
y = [0; 1; 2];
w = [0.1; -0.2; 0.3];

[nll, grad] = neglogli_poiss(w, X, y, []);
expectedNll = sum(exp(X * w) - y .* (X * w));
expectedGrad = full(X' * (exp(X * w) - y));

testCase.verifyEqual(nll, expectedNll, 'AbsTol', 1e-12);
testCase.verifyEqual(grad, expectedGrad, 'AbsTol', 1e-12);
end
