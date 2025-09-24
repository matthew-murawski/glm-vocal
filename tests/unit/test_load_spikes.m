function tests = test_load_spikes
%TEST_LOAD_SPIKES Unit tests for load_spikes.
%   tests = TEST_LOAD_SPIKES() returns function-based tests.

% register local helpers with MATLAB's test harness
tests = functiontests(localfunctions);
end

function testHappyPath(testCase) %#ok<INUSD>
% write a valid spike file for a single neuron
file = writeTempSpikeMat(struct( ...
    'spike_times', (0.1:0.1:0.5)', ...
    'neuron_id', 'n001', ...
    'session_id', 'sessA'));
cleanup = onCleanup(@() deleteIfExists(file));

% load the spikes and confirm all fields round-trip
sp = load_spikes(file);
assert(isequal(sp.neuron_id, 'n001'));
assert(isequal(sp.session_id, 'sessA'));
assert(isequal(sp.spike_times, (0.1:0.1:0.5)'));
clear cleanup
end

function testUnsortedSpikeTimes(testCase)
% create a spike file with unsorted times
file = writeTempSpikeMat(struct( ...
    'spike_times', [0.3; 0.1; 0.2], ...
    'neuron_id', 'n002', ...
    'session_id', 'sessB'));
cleanup = onCleanup(@() deleteIfExists(file));

% expect a warning and verify times come back sorted
verifyWarning(testCase, @() load_spikes(file), 'glm:UnsortedSpikeTimes');
sp = load_spikes(file);
assert(issorted(sp.spike_times));
clear cleanup
end

function testNonFiniteSpikeTimes(testCase)
% insert a nan value to trigger validation failure
file = writeTempSpikeMat(struct( ...
    'spike_times', [0.1; NaN; 0.2], ...
    'neuron_id', 'n003', ...
    'session_id', 'sessC'));
cleanup = onCleanup(@() deleteIfExists(file));

% loader should abort when non-finite spikes appear
verifyError(testCase, @() load_spikes(file), 'glm:InvalidSpikeTimes');
clear cleanup
end

function testMissingFields(testCase)
% omit required metadata to ensure the guard fires
file = writeTempSpikeMat(struct( ...
    'spike_times', (0.1:0.1:0.3)'));
cleanup = onCleanup(@() deleteIfExists(file));

% missing neuron/session fields should raise an error
verifyError(testCase, @() load_spikes(file), 'glm:InvalidSpikesStruct');
clear cleanup
end

function file = writeTempSpikeMat(data)
% helper to persist a struct as a temp mat file
file = [tempname, '.mat'];
save(file, '-struct', 'data');
end

function deleteIfExists(file)
% clean up temp files created during the test run
if exist(file, 'file') == 2
    delete(file);
end
end
