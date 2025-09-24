function tests = test_validate_inputs
%TEST_VALIDATE_INPUTS Unit tests for validate_inputs gatekeeper.
%   tests = TEST_VALIDATE_INPUTS() returns function-based tests.

% register local test functions with MATLAB's test harness
tests = functiontests(localfunctions);
end

function testValidInputs(testCase) %#ok<INUSD>
% craft aligned spike and event data within the allowed window
sp = struct('spike_times', [0.1; 0.2; 0.3]);
ev = struct('t_on', {0.05, 0.25}, 't_off', {0.06, 0.3});

% should run without throwing when all conditions satisfied
validate_inputs(sp, ev, 0.5, 0.1);
end

function testInvalidSpikeOrdering(testCase)
% create spikes that violate monotonic requirement
sp = struct('spike_times', [0.2; 0.1]);
ev = struct('t_on', {0.0}, 't_off', {0.05});

% expect a spike-specific error when sequence decreases
verifyError(testCase, @() validate_inputs(sp, ev, 1.0, 0.1), 'glm:InvalidSpikeTimes');
end

function testInvalidSpikeFinite(testCase)
% inject a nan spike to trigger the finite check
sp = struct('spike_times', [0.1; NaN]);
ev = struct('t_on', {0.0}, 't_off', {0.05});

% loader should flag nan spikes as invalid
verifyError(testCase, @() validate_inputs(sp, ev, 0.5, 0.1), 'glm:InvalidSpikeTimes');
end

function testInvalidEventOrdering(testCase)
% build an event with inverted on/off times
sp = struct('spike_times', [0.1; 0.2]);
ev = struct('t_on', {0.2}, 't_off', {0.1});

% validation must fail when t_off precedes t_on
verifyError(testCase, @() validate_inputs(sp, ev, 0.5, 0.1), 'glm:InvalidLabels');
end

function testEventTimesOutOfBounds(testCase)
% push events beyond the allowed session window after buffer
sp = struct('spike_times', [0.1; 0.2]);
ev = struct('t_on', {-0.2}, 't_off', {-0.1});

% with a small buffer the negative event should trip the guard
verifyError(testCase, @() validate_inputs(sp, ev, 1.0, 0.05), 'glm:InvalidLabels');
end

function testSpikeTimesOutOfBounds(testCase)
% place a spike past the permissible session horizon
sp = struct('spike_times', [0.1; 1.2]);
ev = struct('t_on', {0.2}, 't_off', {0.25});

% expect spike-specific failure when exceeding session bounds
verifyError(testCase, @() validate_inputs(sp, ev, 1.0, 0.1), 'glm:InvalidSpikeTimes');
end

function testEventsMissingFields(testCase)
% pass a malformed event struct lacking required columns
sp = struct('spike_times', [0.1; 0.2]);
ev = struct('t_on', {0.1});

% validator should highlight the schema problem
verifyError(testCase, @() validate_inputs(sp, ev, 1.0, 0.1), 'glm:InvalidLabels');
end
