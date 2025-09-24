function tests = test_load_spikes
%TEST_LOAD_SPIKES Unit tests for load_spikes.
%   tests = TEST_LOAD_SPIKES() returns function-based tests.

tests = functiontests(localfunctions);
end

function testHappyPath(testCase) %#ok<INUSD>
file = writeTempSpikeMat(struct( ...
    'spike_times', (0.1:0.1:0.5)', ...
    'neuron_id', 'n001', ...
    'session_id', 'sessA'));
cleanup = onCleanup(@() deleteIfExists(file));

sp = load_spikes(file);
assert(isequal(sp.neuron_id, 'n001'));
assert(isequal(sp.session_id, 'sessA'));
assert(isequal(sp.spike_times, (0.1:0.1:0.5)'));
clear cleanup
end

function testUnsortedSpikeTimes(testCase)
file = writeTempSpikeMat(struct( ...
    'spike_times', [0.3; 0.1; 0.2], ...
    'neuron_id', 'n002', ...
    'session_id', 'sessB'));
cleanup = onCleanup(@() deleteIfExists(file));

verifyWarning(testCase, @() load_spikes(file), 'glm:UnsortedSpikeTimes');
sp = load_spikes(file);
assert(issorted(sp.spike_times));
clear cleanup
end

function testNonFiniteSpikeTimes(testCase)
file = writeTempSpikeMat(struct( ...
    'spike_times', [0.1; NaN; 0.2], ...
    'neuron_id', 'n003', ...
    'session_id', 'sessC'));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_spikes(file), 'glm:InvalidSpikeTimes');
clear cleanup
end

function testMissingFields(testCase)
file = writeTempSpikeMat(struct( ...
    'spike_times', (0.1:0.1:0.3)'));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_spikes(file), 'glm:InvalidSpikesStruct');
clear cleanup
end

function file = writeTempSpikeMat(data)
file = [tempname, '.mat'];
save(file, '-struct', 'data');
end

function deleteIfExists(file)
if exist(file, 'file') == 2
    delete(file);
end
end
