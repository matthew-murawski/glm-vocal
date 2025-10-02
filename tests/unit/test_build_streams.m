function tests = test_build_streams
% exercise build_streams to confirm interval rasterization rules.

tests = functiontests(localfunctions);
end

function stim = makeStim()
t = (0:0.1:0.5)';
stim = struct('t', t, 'dt', 0.1, 'mask', struct('good', true(numel(t), 1)));
end

function testSingleInterval(testCase)
stim = makeStim();
ev = struct('kind', 'perceived', 't_on', 0.05, 't_off', 0.25, 'label', "");

streams = build_streams(ev, stim);

expectedHeard = logical([1; 0; 0; 0; 0; 0]);
expectedProduced = false(numel(stim.t), 1);

testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_spontaneous, expectedProduced);
testCase.verifyEqual(streams.produced_after_heard, expectedProduced);
testCase.verifyEqual(streams.produced_after_produced, expectedProduced);
testCase.verifyEqual(streams.produced_any, expectedProduced);
testCase.verifyEqual(streams.produced_fields, {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'});
end

function testOverlappingIntervals(testCase)
stim = makeStim();
ev = struct( ...
    'kind', {'perceived', 'produced', 'produced', 'perceived', 'produced'}, ...
    't_on', {0.05, 0.10, 0.20, 0.30, 0.40}, ...
    't_off', {0.15, 0.11, 0.25, 0.35, 0.45}, ...
    'label', {"", "", "", "", ""} ...
);

streams = build_streams(ev, stim);

expectedHeard = logical([1; 0; 0; 1; 0; 0]);
expectedAfterHeard = logical([0; 1; 0; 0; 1; 0]);
expectedAfterProduced = logical([0; 0; 1; 0; 0; 0]);
expectedSpont = false(numel(stim.t), 1);

testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_after_heard, expectedAfterHeard);
testCase.verifyEqual(streams.produced_after_produced, expectedAfterProduced);
testCase.verifyEqual(streams.produced_spontaneous, expectedSpont);
testCase.verifyEqual(streams.produced_any, expectedAfterHeard | expectedAfterProduced);
testCase.verifyEqual(streams.produced_fields, {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'});
end

function testEdgeInclusionConventions(testCase)
stim = makeStim();
ev = struct( ...
    'kind', {'perceived', 'produced', 'produced'}, ...
    't_on', {0.10, 0.50, 0.60}, ...
    't_off', {0.20, 0.60, 0.70}, ...
    'label', {"", "", ""} ...
);
ev(4) = struct('kind', 'perceived', 't_on', 0.20, 't_off', 0.20, 'label', "");

streams = build_streams(ev, stim);

expectedHeard = logical([0; 1; 0; 0; 0; 0]);
expectedProduced = logical([0; 0; 0; 0; 0; 1]);
testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_after_heard, expectedProduced);
testCase.verifyEqual(streams.produced_spontaneous, false(numel(stim.t), 1));
testCase.verifyEqual(streams.produced_after_produced, false(numel(stim.t), 1));
testCase.verifyEqual(streams.produced_any, expectedProduced);
testCase.verifyEqual(streams.produced_fields, {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'});
end

function testProducedContextLookback(testCase)
dt = 0.5;
t = (0:dt:10)';
stim = struct('t', t, 'dt', dt, 'mask', struct('good', true(numel(t), 1)));

ev = struct( ...
    'kind', {'perceived', 'produced', 'produced', 'produced'}, ...
    't_on', {1.0, 2.0, 2.5, 8.0}, ...
    't_off', {1.2, 2.1, 2.6, 8.2}, ...
    'label', {"", "", "", ""} ...
);

streams = build_streams(ev, stim);

expectedSpont = false(numel(t), 1);
expectedSpont(17) = true;
expectedAfterHeard = false(numel(t), 1);
expectedAfterHeard(5) = true;
expectedAfterProduced = false(numel(t), 1);
expectedAfterProduced(6) = true;

testCase.verifyEqual(streams.produced_spontaneous, expectedSpont);
testCase.verifyEqual(streams.produced_after_heard, expectedAfterHeard);
testCase.verifyEqual(streams.produced_after_produced, expectedAfterProduced);
testCase.verifyEqual(streams.produced_any, expectedSpont | expectedAfterHeard | expectedAfterProduced);
testCase.verifyEqual(streams.produced_fields, {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'});
end

function testCallTypeMode(testCase)
stim = makeStim();
ev = struct( ...
    'kind', {'produced', 'produced', 'produced', 'produced'}, ...
    't_on', {0.05, 0.15, 0.25, 0.35}, ...
    't_off', {0.06, 0.16, 0.26, 0.36}, ...
    'label', {"phee", "twitter", "trill", "trill-phee"} ...
);

cfg = struct('produced_split_mode', 'call_type');
streams = build_streams(ev, stim, cfg);

expectedFields = {'produced_phee', 'produced_twitter', 'produced_trill', 'produced_trillphee'};
testCase.verifyEqual(streams.produced_fields, expectedFields);
for ii = 1:numel(expectedFields)
    vec = streams.(expectedFields{ii});
    testCase.verifyEqual(sum(vec), 1);
end
end
