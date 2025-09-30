function tests = test_assemble_design_matrix
% section registration
% register the local tests with the matlab runner.
tests = functiontests(localfunctions);
end

function stim = makeStim(nT, dt, goodMask)
% section stimulus helper
% craft a simple stimulus struct with optional good-mask entries.
stim = struct();
stim.dt = dt;
stim.t = (0:dt:dt * (nT - 1))';
if nargin < 3
    goodMask = true(nT, 1);
end
stim.mask = struct('good', logical(goodMask(:)));
end

function cfg = makeCfg()
% section config helper
% provide compact windows suited to the unit tests.
cfg = struct();
cfg.heard_window_s = [0.0, 0.2];
cfg.produced_window_s = [-0.1, 0.2];
cfg.produced_basis = struct('kind', 'raised_cosine', 'n_basis', 3, 'overlap', 1.5, 'normalize', 'l1');
cfg.history_window_s = [0.1, 0.3];
end

function testColumnOrderAndCounts(testCase)
% section column validation
% we verify the assembler creates the expected number of columns and maintains block ordering.
stim = makeStim(6, 0.1);
streams = struct();
streams.heard_any = [0; 1; 0; 1; 0; 0];
streams.produced_any = [1; 0; 0; 1; 0; 1];

states = struct();
states.convo = [1; 0; 1; 0; 1; 0];
states.spon = double(~states.convo);

sps = (0:5)';
cfg = makeCfg();

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);

testCase.verifyTrue(issparse(Xd.X));
testCase.verifyEqual(size(Xd.X, 1), numel(stim.t));
testCase.verifyEqual(size(Xd.X, 2), 12);
testCase.verifyEqual(Xd.y, sps);

colmap = Xd.colmap;
testCase.verifyEqual(colmap.intercept.cols, 1);
testCase.verifyEqual(colmap.heard_any.cols, (2:4));
testCase.verifyEqual(colmap.produced_any.cols, (5:7));
testCase.verifyEqual(colmap.states.cols, (8:9));
testCase.verifyEqual(colmap.states.convo, 8);
testCase.verifyEqual(colmap.states.spon, 9);
testCase.verifyEqual(colmap.spike_history.cols, (10:12));

Xfull = full(Xd.X);
testCase.verifyEqual(Xfull(:, colmap.intercept.cols), ones(6, 1));
testCase.verifyEqual(colmap.produced_any.info.mode, 'raised_cosine');
testCase.verifyEqual(colmap.produced_any.info.basis.n_basis, 3);
testCase.verifyEqual(size(colmap.produced_any.info.basis.matrix, 2), 3);
testCase.verifyEqual(colmap.heard_any.info.mode, 'causal');
testCase.verifyEqual(colmap.spike_history.info.mode, 'history');
end

function testMaskDropsBadRows(testCase)
% section mask behaviour
% we ensure the assembler respects mask.good by removing bad rows from X and y.
goodMask = [true; false; true; false; true];
stim = makeStim(numel(goodMask), 0.1, goodMask);
streams = struct('heard_any', double((1:5)' > 2), 'produced_any', double(mod((1:5)', 2)));
states = struct('convo', double([1; 0; 1; 0; 1]), 'spon', double([0; 1; 0; 1; 0]));
sps = (5:9)';
cfg = makeCfg();

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);

expectedRows = sum(goodMask);
testCase.verifyEqual(size(Xd.X, 1), expectedRows);
testCase.verifyEqual(numel(Xd.y), expectedRows);
testCase.verifyEqual(Xd.y, sps(goodMask));

Xfull = full(Xd.X);
colmap = Xd.colmap;
testCase.verifyEqual(Xfull(:, colmap.intercept.cols), ones(expectedRows, 1));
testCase.verifyEqual(Xfull(:, colmap.states.convo), states.convo(goodMask));
testCase.verifyEqual(Xfull(:, colmap.states.spon), states.spon(goodMask));
end
