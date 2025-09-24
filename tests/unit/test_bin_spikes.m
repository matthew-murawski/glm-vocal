function tests = test_bin_spikes
%TEST_BIN_SPIKES Unit tests for bin_spikes helper.
%   tests = TEST_BIN_SPIKES() returns function-based tests.

% register local test functions with the MATLAB harness
tests = functiontests(localfunctions);
end

function stim = makeStim()
% helper to create a simple stimulus timeline for the tests
stim = struct('t', (0:0.1:0.5)', 'dt', 0.1, 'mask', struct('good', true(6, 1)));
end

function testCountsMatchSpikeTotal(testCase) %#ok<INUSD>
% prepare spikes across several bins
spikes = [0.05; 0.12; 0.31];
stim = makeStim();

% bin them and confirm the total count matches
sps = bin_spikes(spikes, stim);
assert(sum(sps) == numel(spikes));
end

function testEdgeSpikes(testCase)
% place spikes exactly on bin edges to validate histogram behaviour
spikes = [0.0; 0.1; 0.2; 0.5];
stim = makeStim();

% counts should reflect left-inclusive edges, right-exclusive final
sps = bin_spikes(spikes, stim);
expected = [1; 1; 1; 0; 0; 1];
assert(isequal(sps, expected));
end

function testEmptySpikes(testCase) %#ok<INUSD>
% running on empty input should yield zeros for each bin
spikes = [];
stim = makeStim();

% verify zero counts without errors
sps = bin_spikes(spikes, stim);
assert(isequal(sps, zeros(numel(stim.t), 1)));
end
