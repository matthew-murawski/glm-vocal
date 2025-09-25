function tests = test_predict_and_unpack
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testPredictRateClipping(testCase)
% section rate prediction
% verify mu = exp(X*w) with clipping to avoid overflow.
X = [1 1; -1 -1; 0 0];
w = [40; 40];
mu = predict_rate(X, w);

expected = exp([50; -50; 0]);
testCase.verifyEqual(mu, expected, 'AbsTol', 1e-12);
end

function testUnpackParamsMapping(testCase)
% section parameter unpacking
% ensure weights are sliced into named kernels according to the column map metadata.
w = (1:8)' * 0.1;
wmap = struct('w', w);

colmap = struct();
colmap.intercept = struct('cols', 1, 'name', 'intercept');
colmap.heard_any = struct('cols', 2:3, 'info', struct('lag_bins', (0:1)', 'lag_times_s', (0:1)' * 0.01, 'window_s', [0, 0.01], 'mode', 'causal'));
colmap.produced_any = struct('cols', 4:5, 'info', struct('lag_bins', (-1:0)', 'lag_times_s', (-1:0)' * 0.01, 'window_s', [-0.01, 0], 'mode', 'symmetric'));
colmap.spike_history = struct('cols', 6, 'info', struct('lag_bins', 1, 'lag_times_s', 0.01, 'window_s', [0.01, 0.01], 'mode', 'history'));
colmap.states = struct('cols', [7 8], 'names', {{'convo', 'spon'}}, 'convo', 7, 'spon', 8);

kernels = unpack_params(wmap, colmap, struct(), struct());

testCase.verifyEqual(kernels.intercept, w(1));

expectedHeard = w(2:3);
testCase.verifyEqual(kernels.heard_any.weights, expectedHeard);
testCase.verifyEqual(kernels.heard_any.mode, 'causal');
testCase.verifyEqual(kernels.heard_any.lag_bins, (0:1)');

testCase.verifyEqual(kernels.produced_any.weights, w(4:5));
testCase.verifyEqual(kernels.produced_any.mode, 'symmetric');

testCase.verifyEqual(kernels.spike_history.weights, w(6));

testCase.verifyEqual(kernels.states.weights, w(7:8));
testCase.verifyEqual(kernels.states.convo, w(7));
testCase.verifyEqual(kernels.states.spon, w(8));
end
