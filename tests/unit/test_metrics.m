function tests = test_metrics
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testMatchesKnownValues(testCase)
% section known case
% verify metrics against manually computed values on a small dataset.
y = [0; 1; 2; 3];
mu = [0.5; 1.2; 1.8; 2.5];

out = metrics(y, mu);

% manual calculations
LL_const = gammaln(y + 1);
nllExpected = sum(mu - y .* log(mu) + LL_const);

baseline = mean(y);
muNull = baseline * ones(size(y));
LLModel = sum(y .* log(mu) - mu);
LLNull = sum(y .* log(muNull) - muNull);
mask = y > 0;
LLSat = sum(y(mask) .* log(y(mask)) - y(mask));

pseudoExpected = 1 - (LLSat - LLModel) / (LLSat - LLNull);

devModel = 2 * (LLSat - LLModel);
devNull = 2 * (LLSat - LLNull);
DevExpected = 1 - devModel / devNull;

compareStructs(testCase, out, struct('nll', nllExpected, 'pseudoR2', pseudoExpected, 'deviance_explained', DevExpected));
end

function testHandlesZeroVariance(testCase)
% section degenerate baseline
% ensure metrics gracefully handle cases where the null deviance is zero.
y = [1; 1; 1];
mu = [1; 1; 1];

out = metrics(y, mu);

testCase.verifyTrue(isnan(out.pseudoR2));
testCase.verifyTrue(isnan(out.deviance_explained));
end

function compareStructs(testCase, actual, expected)
fields = fieldnames(expected);
for ii = 1:numel(fields)
    field = fields{ii};
    expVal = expected.(field);
    actVal = actual.(field);
    testCase.verifyEqual(actVal, expVal, 'AbsTol', 1e-12);
end
end
