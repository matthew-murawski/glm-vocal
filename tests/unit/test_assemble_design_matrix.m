function tests = test_assemble_design_matrix
% section registration
% register the local tests with the matlab runner.
tests = functiontests(localfunctions);
end

function stim = makeStim(nT, dt, goodMask)
stim = struct();
stim.dt = dt;
stim.t = (0:dt:dt * (nT - 1))';
if nargin < 3
    goodMask = true(nT, 1);
end
stim.mask = struct('good', logical(goodMask(:)));
end

function cfg = makeCfg()
cfg = struct();
cfg.heard_window_s = [0.0, 0.2];
cfg.heard_basis = struct('kind', 'raised_cosine', 'n_basis', 3, 'overlap', 1.5, 'normalize', 'l1');
cfg.produced_window_s = [-0.1, 0.2];
cfg.produced_basis = struct('kind', 'raised_cosine', 'n_basis', 3, 'overlap', 1.5, 'normalize', 'l1');
cfg.history_window_s = [0.1, 0.3];
end

function testColumnOrderAndCounts(testCase)
stim = makeStim(6, 0.1);
streams = struct();
streams.heard_any = [0; 1; 0; 1; 0; 0];
streams.heard_addressed = streams.heard_any;
streams.heard_overheard = zeros(size(streams.heard_any));
streams.heard_fields = {'heard_addressed', 'heard_overheard'};
streams.produced_spontaneous = [1; 0; 0; 0; 1; 0];
streams.produced_after_heard = [0; 1; 0; 0; 0; 1];
streams.produced_after_produced = [0; 0; 1; 0; 0; 0];
streams.produced_fields = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
streams.produced_any = streams.produced_spontaneous | streams.produced_after_heard | streams.produced_after_produced;

states = struct();
states.convo = [1; 0; 1; 0; 1; 0];
states.spon = double(~states.convo);

sps = (0:5)';
cfg = makeCfg();

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);

testCase.verifyTrue(issparse(Xd.X));
testCase.verifyEqual(size(Xd.X, 1), numel(stim.t));
testCase.verifyEqual(size(Xd.X, 2), 21);
testCase.verifyEqual(Xd.y, sps);

colmap = Xd.colmap;
testCase.verifyEqual(colmap.heard_fields, streams.heard_fields);
testCase.verifyEqual(colmap.produced_fields, streams.produced_fields);
testCase.verifyEqual(colmap.intercept.cols, 1);
testCase.verifyEqual(colmap.heard_addressed.cols, (2:4));
testCase.verifyEqual(colmap.heard_overheard.cols, (5:7));
testCase.verifyEqual(colmap.produced_spontaneous.cols, (8:10));
testCase.verifyEqual(colmap.produced_after_heard.cols, (11:13));
testCase.verifyEqual(colmap.produced_after_produced.cols, (14:16));
testCase.verifyEqual(colmap.states.cols, (17:18));
testCase.verifyEqual(colmap.states.convo, 17);
testCase.verifyEqual(colmap.states.spon, 18);
testCase.verifyEqual(colmap.spike_history.cols, (19:21));
end

function testMaskDropsBadRows(testCase)
goodMask = [true; false; true; false; true];
stim = makeStim(numel(goodMask), 0.1, goodMask);
streams = struct();
streams.heard_any = double((1:5)' > 2);
streams.heard_addressed = streams.heard_any;
streams.heard_overheard = zeros(size(streams.heard_any));
streams.heard_fields = {'heard_addressed', 'heard_overheard'};
streams.produced_spontaneous = double(mod((1:5)', 2) == 0);
streams.produced_after_heard = double(mod((1:5)', 2) == 1);
streams.produced_after_produced = zeros(5, 1);
streams.produced_fields = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
streams.produced_any = streams.produced_spontaneous | streams.produced_after_heard | streams.produced_after_produced;
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

function testExcludePredictorsSkipsBlocks(testCase)
stim = makeStim(6, 0.1);
streams = struct();
streams.heard_any = [0; 1; 0; 1; 0; 0];
streams.heard_addressed = streams.heard_any;
streams.heard_overheard = zeros(size(streams.heard_any));
streams.heard_fields = {'heard_addressed', 'heard_overheard'};
streams.produced_spontaneous = [1; 0; 0; 0; 1; 0];
streams.produced_after_heard = [0; 1; 0; 0; 0; 1];
streams.produced_after_produced = [0; 0; 1; 0; 0; 0];
streams.produced_fields = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
streams.produced_any = streams.produced_spontaneous | streams.produced_after_heard | streams.produced_after_produced;

states = struct();
states.convo = [1; 0; 1; 0; 1; 0];
states.spon = double(~states.convo);

sps = (0:5)';
cfg = makeCfg();
cfg.exclude_predictors = {'states', 'spike_history'};

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);

testCase.verifyEqual(size(Xd.X, 2), 16);
testCase.verifyFalse(isfield(Xd.colmap, 'states'));
testCase.verifyFalse(isfield(Xd.colmap, 'spike_history'));
testCase.verifyTrue(isfield(Xd.colmap, 'heard_addressed'));
testCase.verifyTrue(isfield(Xd.colmap, 'heard_overheard'));
testCase.verifyEqual(Xd.y, sps);
end

function testCallTypeColumns(testCase)
stim = makeStim(4, 0.1);
streams = struct();
streams.heard_any = zeros(4, 1);
streams.heard_addressed = zeros(4, 1);
streams.heard_overheard = zeros(4, 1);
streams.heard_fields = {'heard_addressed', 'heard_overheard'};
streams.produced_phee = [1; 0; 0; 0];
streams.produced_twitter = [0; 1; 0; 0];
streams.produced_trill = [0; 0; 1; 0];
streams.produced_trillphee = [0; 0; 0; 1];
streams.produced_fields = {'produced_phee', 'produced_twitter', 'produced_trill', 'produced_trillphee'};
streams.produced_any = streams.produced_phee | streams.produced_twitter | streams.produced_trill | streams.produced_trillphee;

states = struct('convo', zeros(4, 1), 'spon', ones(4, 1));
sps = (1:4)';
cfg = makeCfg();

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);

colmap = Xd.colmap;
testCase.verifyTrue(isfield(colmap, 'produced_phee'));
testCase.verifyEqual(colmap.produced_fields, streams.produced_fields);
end
