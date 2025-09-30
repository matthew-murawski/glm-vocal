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
testCase.verifyEqual(streams.produced_any, expectedProduced);
end

function testOverlappingIntervals(testCase)
% overlapping events of each kind should union correctly across bins.

%% craft perceived and produced events with overlapping spans
stim = makeStim();
perceivedEv = struct( ...
    'kind', {'perceived', 'perceived'}, ...
    't_on', {0.05, 0.15}, ...
    't_off', {0.25, 0.35}, ...
    'label', {"", ""} ...
);
producedEv = struct( ...
    'kind', {'produced', 'produced'}, ...
    't_on', {0.10, 0.30}, ...
    't_off', {0.40, 0.50}, ...
    'label', {"", ""} ...
);
ev = [perceivedEv(:); producedEv(:)];

%% run the stream builder
streams = build_streams(ev, stim);

%% confirm the perceived and produced impulse streams align with expectations
testCase.verifyEqual(streams.heard_any, logical([1; 1; 0; 0; 0; 0]));
testCase.verifyEqual(streams.produced_any, logical([0; 1; 0; 1; 0; 0]));
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
testCase.verifyEqual(streams.produced_any, expectedProduced);
end
