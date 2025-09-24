function tests = test_load_labels
%TEST_LOAD_LABELS Unit tests for load_labels.
%   tests = TEST_LOAD_LABELS() returns function-based tests.

% expose local test helpers to the MATLAB harness
tests = functiontests(localfunctions);
end

function testLoadMatHappy(testCase)
% assemble a proper events struct with two entries
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced', 'perceived'}, ...
        'onset', {0.1, 0.5}, ...
        'offset', {0.2, 0.6}, ...
        'label', {'callA', 'callB'})));
cleanup = onCleanup(@() deleteIfExists(file));

% run the loader and make sure data round-trips cleanly
labels = load_labels(file);
assert(numel(labels) == 2);
assert(strcmp(labels(1).kind, 'produced'));
assert(strcmp(labels(2).kind, 'perceived'));
assert(labels(1).t_on == 0.1 && labels(1).t_off == 0.2);
assert(labels(2).t_on == 0.5 && labels(2).t_off == 0.6);
assert(labels(1).label == "callA");
assert(labels(2).label == "callB");
clear cleanup
end

function testUnsupportedExtension(testCase)
% create a bogus txt file to check the extension guard
file = [tempname, '.txt'];
fid = fopen(file, 'w');
assert(fid ~= -1);
fprintf(fid, 'dummy');
fclose(fid);
cleanup = onCleanup(@() deleteIfExists(file));

% loader should reject unsupported formats immediately
verifyError(testCase, @() load_labels(file), 'glm:UnsupportedLabelFormat');
clear cleanup
end

function testMissingField(testCase)
% drop the offset field to trigger schema validation
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced'}, ...
        'onset', {0.1}, ...
        'label', {'x'})));
cleanup = onCleanup(@() deleteIfExists(file));

% expect an invalid-struct error when required fields are absent
verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelsStruct');
clear cleanup
end

function testInvalidKind(testCase)
% provide an unsupported kind label to exercise the guard
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'other'}, ...
        'onset', {0.0}, ...
        'offset', {0.1}, ...
        'label', {'unknown'})));
cleanup = onCleanup(@() deleteIfExists(file));

% loader should flag invalid kinds before returning data
verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelKind');
clear cleanup
end

function testInvalidTimes(testCase)
% capture a negative duration event to test time validation
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced'}, ...
        'onset', {0.2}, ...
        'offset', {0.1}, ...
        'label', {'bad'})));
cleanup = onCleanup(@() deleteIfExists(file));

% offset preceding onset should raise an invalid-times error
verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelTimes');
clear cleanup
end

function file = writeTempMat(data)
% helper to stash an events struct as a temporary mat file
file = [tempname, '.mat'];
save(file, '-struct', 'data');
end

function deleteIfExists(file)
% clean up any temporary files emitted during testing
if exist(file, 'file') == 2
    delete(file);
end
end
