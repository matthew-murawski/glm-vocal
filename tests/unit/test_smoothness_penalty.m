function tests = test_smoothness_penalty
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function colmap = makeColmap()
% section helper
% craft a representative column map to drive the penalty assembly tests.
colmap = struct();
colmap.intercept = struct('cols', 1, 'name', 'intercept');
colmap.heard_any = struct('cols', 2:5, 'info', struct());
colmap.produced_any = struct('cols', 6:9, 'info', struct());
colmap.states = struct('cols', 10:11, 'names', {{'convo', 'spon'}}, 'convo', 10, 'spon', 11);
colmap.spike_history = struct('cols', 12:15, 'info', struct());
end

function testSecondDifferenceStencils(testCase)
% section stencil validation
% confirm that second-difference rows are generated for each kernel block while leaving scalar columns unpenalised.
colmap = makeColmap();
[D, Dmap] = smoothness_penalty(colmap, struct());

testCase.verifyEqual(size(D), [6, 15]);
heardRows = Dmap.heard_any.rows;
expectedHeard = [1 -2 1 0; 0 1 -2 1];
fullHeard = full(D(heardRows, colmap.heard_any.cols));
testCase.verifyEqual(fullHeard, expectedHeard);
testCase.verifyEqual(full(D(:, colmap.intercept.cols)), zeros(6, 1));
testCase.verifyEqual(full(D(:, colmap.states.cols)), zeros(6, 2));
testCase.verifyEqual(full(D(Dmap.produced_any.rows, colmap.produced_any.cols)), expectedHeard);
testCase.verifyEqual(full(D(Dmap.spike_history.rows, colmap.spike_history.cols)), expectedHeard);
end

function testConfigDisablesBlocks(testCase)
% section config behaviour
% ensure the configuration can disable penalties for individual blocks.
colmap = makeColmap();
cfg = struct('produced_any', false);
[D, Dmap] = smoothness_penalty(colmap, cfg);
testCase.verifyFalse(isfield(Dmap, 'produced_any'));
testCase.verifyEqual(size(D), [4, 15]);
testCase.verifyEqual(full(D(:, colmap.produced_any.cols)), zeros(4, numel(colmap.produced_any.cols)));
end

function testLambdaScalingPerBlock(testCase)
% section lambda scaling
% verify per-block lambdas scale second-difference stencils as expected while allowing zero to disable a block.
colmap = makeColmap();
lambdaCfg = struct('heard', 4, 'produced', 9, 'history', 0);
[D, Dmap] = smoothness_penalty(colmap, struct(), lambdaCfg);

sqrtHeard = sqrt(lambdaCfg.heard);
sqrtProduced = sqrt(lambdaCfg.produced);

heardBlock = full(D(Dmap.heard_any.rows, colmap.heard_any.cols)) ./ sqrtHeard;
expectedHeard = [1 -2 1 0; 0 1 -2 1];
testCase.verifyEqual(heardBlock, expectedHeard);

producedBlock = full(D(Dmap.produced_any.rows, colmap.produced_any.cols)) ./ sqrtProduced;
testCase.verifyEqual(producedBlock, expectedHeard);

testCase.verifyTrue(isfield(Dmap, 'spike_history'));
testCase.verifyEqual(numel(Dmap.spike_history.rows), 0);
testCase.verifyEqual(size(D), [4, 15]);
end

function testScalarLambdaAppliesUniformScaling(testCase)
% section scalar lambda
% confirm a scalar lambda scales every penalised block uniformly.
colmap = makeColmap();
lambda = 0.25;
[D, Dmap] = smoothness_penalty(colmap, struct(), lambda);

scale = sqrt(lambda);
heardBlock = full(D(Dmap.heard_any.rows, colmap.heard_any.cols)) ./ scale;
expectedHeard = [1 -2 1 0; 0 1 -2 1];
testCase.verifyEqual(heardBlock, expectedHeard);
producedBlock = full(D(Dmap.produced_any.rows, colmap.produced_any.cols)) ./ scale;
testCase.verifyEqual(producedBlock, expectedHeard);
historyBlock = full(D(Dmap.spike_history.rows, colmap.spike_history.cols)) ./ scale;
testCase.verifyEqual(historyBlock, expectedHeard);
end
