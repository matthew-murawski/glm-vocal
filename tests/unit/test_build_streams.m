function tests = test_build_streams
% exercise build_streams to confirm interval rasterization rules.

%% register the test cases
% we expose the local functions to matlab's function-based test harness.
tests = functiontests(localfunctions);
end

function stim = makeStim()
% build a reusable stimulus timeline for the tests.

%% define the timeline grid
% we keep the grid simple so expected logical vectors stay readable.
t = (0:0.1:0.5)';
stim = struct('t', t, 'dt', 0.1, 'mask', struct('good', true(numel(t), 1)));
end

function testSingleInterval(testCase)
% single perceived interval should mark the overlapping bins only.

%% set up a perceived interval spanning multiple bins
% the event begins mid-bin and ends midway through the third bin.
stim = makeStim();
ev = struct('kind', 'perceived', 't_on', 0.05, 't_off', 0.25, 'label', "");

%% run the stream builder
streams = build_streams(ev, stim);

%% confirm the impulse aligns with the onset bin only
expectedHeard = logical([1; 0; 0; 0; 0; 0]);
expectedProduced = false(numel(stim.t), 1);
testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_spontaneous, expectedProduced);
testCase.verifyEqual(streams.produced_after_heard, expectedProduced);
testCase.verifyEqual(streams.produced_after_produced, expectedProduced);
testCase.verifyEqual(streams.produced_any, expectedProduced);
end

function testOverlappingIntervals(testCase)
% overlapping events of each kind should union correctly across bins.

%% craft perceived and produced events with overlapping spans
stim = makeStim();
ev = struct( ...
    'kind', {'perceived', 'produced', 'produced', 'perceived', 'produced'}, ...
    't_on', {0.05, 0.10, 0.20, 0.30, 0.40}, ...
    't_off', {0.15, 0.11, 0.25, 0.35, 0.45}, ...
    'label', {"", "", "", "", ""} ...
);

%% run the stream builder
streams = build_streams(ev, stim);

%% confirm the perceived and produced impulse streams align with expectations
expectedHeard = logical([1; 0; 0; 1; 0; 0]);
expectedAfterHeard = logical([0; 1; 0; 0; 1; 0]);
expectedAfterProduced = logical([0; 0; 1; 0; 0; 0]);
expectedSpont = false(numel(stim.t), 1);

testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_after_heard, expectedAfterHeard);
testCase.verifyEqual(streams.produced_after_produced, expectedAfterProduced);
testCase.verifyEqual(streams.produced_spontaneous, expectedSpont);
testCase.verifyEqual(streams.produced_any, expectedAfterHeard | expectedAfterProduced);
end

function testEdgeInclusionConventions(testCase)
% check inclusive start, exclusive end, and far-edge handling.

%% assemble events probing the bin boundaries
stim = makeStim();
ev = struct( ...
    'kind', {'perceived', 'produced', 'produced'}, ...
    't_on', {0.10, 0.50, 0.60}, ...
    't_off', {0.20, 0.60, 0.70}, ...
    'label', {"", "", ""} ...
);
ev(4) = struct('kind', 'perceived', 't_on', 0.20, 't_off', 0.20, 'label', "");

%% run the stream builder
streams = build_streams(ev, stim);

%% validate inclusive/exclusive behaviour and out-of-range handling
expectedHeard = logical([0; 1; 0; 0; 0; 0]);
expectedProduced = logical([0; 0; 0; 0; 0; 1]);
testCase.verifyEqual(streams.heard_any, expectedHeard);
testCase.verifyEqual(streams.produced_after_heard, expectedProduced);
testCase.verifyEqual(streams.produced_spontaneous, false(numel(stim.t), 1));
testCase.verifyEqual(streams.produced_after_produced, false(numel(stim.t), 1));
testCase.verifyEqual(streams.produced_any, expectedProduced);
end

function testProducedContextLookback(testCase)
% produced calls should respect the five second lookback window when categorised.

%% build a longer timeline with distinct produced contexts
dt = 0.5;
t = (0:dt:10)';
stim = struct('t', t, 'dt', dt, 'mask', struct('good', true(numel(t), 1)));

ev = struct( ...
    'kind', {'perceived', 'produced', 'produced', 'produced'}, ...
    't_on', {1.0, 2.0, 2.5, 8.0}, ...
    't_off', {1.2, 2.1, 2.6, 8.2}, ...
    'label', {"", "", "", ""} ...
);

%% run the stream builder
streams = build_streams(ev, stim);

%% check each produced context is assigned to the appropriate stream
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
end
