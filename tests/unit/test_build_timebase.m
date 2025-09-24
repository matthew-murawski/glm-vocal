function tests = test_build_timebase
%TEST_BUILD_TIMEBASE Unit tests for the build_timebase helper.
%   tests = TEST_BUILD_TIMEBASE() returns function-based tests.

% register local test functions with MATLAB's unit test harness
tests = functiontests(localfunctions);
end

function testCoversLastObservation(testCase) %#ok<INUSD>
% create spikes and events ending at different times
sp = struct('spike_times', [0.05; 0.1; 0.19]);
ev = struct('t_on', {0.15}, 't_off', {0.22});

% build the timebase and assert it covers the furthest timestamp
stim = build_timebase(ev, sp, 0.02);
assert(stim.t(end) >= 0.22);
end

function testRoundingBehavior(testCase)
% set up a last event slightly beyond an integer multiple of dt
sp = struct('spike_times', [0.0; 0.05]);
ev = struct('t_on', {0.0}, 't_off', {0.101});

% verify rounding logic stays within half-step tolerance below T
stim = build_timebase(ev, sp, 0.05);
expectedBins = floor(0.101 / 0.05 + 0.5) + 1;
assert(numel(stim.t) == expectedBins);
assert(stim.t(end) <= 0.101 + 0.025 + 1e-12);
end
