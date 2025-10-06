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
colmap.heard_fields = {'heard_addressed', 'heard_overheard'};
colmap.heard_addressed = struct('cols', 2:5, 'info', struct());
colmap.heard_overheard = struct('cols', 6:9, 'info', struct());
colmap.produced_spontaneous = struct('cols', 10:13, 'info', struct());
colmap.produced_after_heard = struct('cols', 14:17, 'info', struct());
colmap.produced_after_produced = struct('cols', 18:21, 'info', struct());
colmap.states = struct('cols', 22:23, 'names', {{'convo', 'spon'}}, 'convo', 22, 'spon', 23);
colmap.spike_history = struct('cols', 24:27, 'info', struct());
end

function testSecondDifferenceStencils(testCase)
% section stencil validation
% confirm that second-difference rows are generated for each kernel block while leaving scalar columns unpenalised.
colmap = makeColmap();
[D, Dmap] = smoothness_penalty(colmap, struct());

testCase.verifyEqual(size(D), [12, 27]);
heardRows = Dmap.heard_addressed.rows;
expectedHeard = [1 -2 1 0; 0 1 -2 1];
fullHeard = full(D(heardRows, colmap.heard_addressed.cols));
testCase.verifyEqual(fullHeard, expectedHeard);
overheardRows = Dmap.heard_overheard.rows;
fullOverheard = full(D(overheardRows, colmap.heard_overheard.cols));
testCase.verifyEqual(fullHeard, expectedHeard);
testCase.verifyEqual(fullOverheard, expectedHeard);
testCase.verifyEqual(full(D(:, colmap.intercept.cols)), zeros(12, 1));
testCase.verifyEqual(full(D(:, colmap.states.cols)), zeros(12, 2));
testCase.verifyEqual(full(D(Dmap.produced_spontaneous.rows, colmap.produced_spontaneous.cols)), expectedHeard);
testCase.verifyEqual(full(D(Dmap.produced_after_heard.rows, colmap.produced_after_heard.cols)), expectedHeard);
testCase.verifyEqual(full(D(Dmap.produced_after_produced.rows, colmap.produced_after_produced.cols)), expectedHeard);
testCase.verifyEqual(full(D(Dmap.spike_history.rows, colmap.spike_history.cols)), expectedHeard);
end

function testConfigDisablesBlocks(testCase)
% section config behaviour
% ensure the configuration can disable penalties for individual blocks.
colmap = makeColmap();
cfg = struct('produced_after_heard', false);
[D, Dmap] = smoothness_penalty(colmap, cfg);
testCase.verifyFalse(isfield(Dmap, 'produced_after_heard'));
testCase.verifyEqual(size(D), [10, 27]);
testCase.verifyEqual(full(D(:, colmap.produced_after_heard.cols)), zeros(10, numel(colmap.produced_after_heard.cols)));
end

function testLambdaScalingPerBlock(testCase)
% section lambda scaling
% verify per-block lambdas scale second-difference stencils as expected while allowing zero to disable a block.
colmap = makeColmap();
lambdaCfg = struct('heard', 4, 'produced', 9, 'history', 0);
[D, Dmap] = smoothness_penalty(colmap, struct(), lambdaCfg);

sqrtHeard = sqrt(lambdaCfg.heard);
sqrtProduced = sqrt(lambdaCfg.produced);

expectedHeard = [1 -2 1 0; 0 1 -2 1];
heardBlock = full(D(Dmap.heard_addressed.rows, colmap.heard_addressed.cols)) ./ sqrtHeard;
overheardBlock = full(D(Dmap.heard_overheard.rows, colmap.heard_overheard.cols)) ./ sqrtHeard;
testCase.verifyEqual(heardBlock, expectedHeard);
testCase.verifyEqual(overheardBlock, expectedHeard);

spontBlock = full(D(Dmap.produced_spontaneous.rows, colmap.produced_spontaneous.cols)) ./ sqrtProduced;
afterHeardBlock = full(D(Dmap.produced_after_heard.rows, colmap.produced_after_heard.cols)) ./ sqrtProduced;
afterProducedBlock = full(D(Dmap.produced_after_produced.rows, colmap.produced_after_produced.cols)) ./ sqrtProduced;
testCase.verifyEqual(spontBlock, expectedHeard);
testCase.verifyEqual(afterHeardBlock, expectedHeard);
testCase.verifyEqual(afterProducedBlock, expectedHeard);

testCase.verifyTrue(isfield(Dmap, 'spike_history'));
testCase.verifyEqual(numel(Dmap.spike_history.rows), 0);
testCase.verifyEqual(size(D), [10, 27]);
end

function testScalarLambdaAppliesUniformScaling(testCase)
% section scalar lambda
% confirm a scalar lambda scales every penalised block uniformly.
colmap = makeColmap();
lambda = 0.25;
[D, Dmap] = smoothness_penalty(colmap, struct(), lambda);

scale = sqrt(lambda);
expectedHeard = [1 -2 1 0; 0 1 -2 1];
heardBlock = full(D(Dmap.heard_addressed.rows, colmap.heard_addressed.cols)) ./ scale;
overheardBlock = full(D(Dmap.heard_overheard.rows, colmap.heard_overheard.cols)) ./ scale;
spontBlock = full(D(Dmap.produced_spontaneous.rows, colmap.produced_spontaneous.cols)) ./ scale;
afterHeardBlock = full(D(Dmap.produced_after_heard.rows, colmap.produced_after_heard.cols)) ./ scale;
afterProducedBlock = full(D(Dmap.produced_after_produced.rows, colmap.produced_after_produced.cols)) ./ scale;
historyBlock = full(D(Dmap.spike_history.rows, colmap.spike_history.cols)) ./ scale;

testCase.verifyEqual(heardBlock, expectedHeard);
testCase.verifyEqual(overheardBlock, expectedHeard);
testCase.verifyEqual(spontBlock, expectedHeard);
testCase.verifyEqual(afterHeardBlock, expectedHeard);
testCase.verifyEqual(afterProducedBlock, expectedHeard);
testCase.verifyEqual(historyBlock, expectedHeard);
testCase.verifyEqual(size(D), [12, 27]);
end
