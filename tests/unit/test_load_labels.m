function tests = test_load_labels
% exercise load_labels across mat and txt inputs.

%% register local tests with the matlab harness
% we expose the local test functions to function-based testing.
tests = functiontests(localfunctions);
end

function testLoadMatHappy(testCase)
% assemble a proper events struct with two entries.
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced', 'perceived'}, ...
        'onset', {0.1, 0.5}, ...
        'offset', {0.2, 0.6}, ...
        'label', {'callA', 'callB'})));
cleanup = onCleanup(@() deleteIfExists(file));

labels = load_labels(file);
verifyEqual(testCase, numel(labels), 2);
verifyEqual(testCase, labels(1).kind, 'produced');
verifyEqual(testCase, labels(2).kind, 'perceived');
verifyEqual(testCase, labels(1).t_on, 0.1);
verifyEqual(testCase, labels(1).t_off, 0.2);
verifyEqual(testCase, labels(1).label, "callA");
verifyEqual(testCase, labels(2).label, "callB");
clear cleanup
end

function testLoadTxtWithKinds(testCase)
% parse a txt file whose label column contains produced/perceived kinds.
file = writeTempTxt({ ...
    sprintf('0.00\t0.10\tproduced'), ...
    sprintf('0.20\t0.30\tperceived') ...
});
cleanup = onCleanup(@() deleteIfExists(file));

labels = load_labels(file);
verifyEqual(testCase, numel(labels), 2);
verifyEqual(testCase, labels(1).kind, 'produced');
verifyEqual(testCase, labels(2).kind, 'perceived');
verifyEqual(testCase, labels(1).label, "produced");
verifyEqual(testCase, labels(2).label, "perceived");
clear cleanup
end

function testLoadTxtWithDefaultKind(testCase)
% ensure txt files without a label column can use a default kind override.
file = writeTempTxt({ ...
    sprintf('0.05\t0.15'), ...
    sprintf('0.20\t0.40') ...
});
cleanup = onCleanup(@() deleteIfExists(file));

labels = load_labels(file, 'produced');
verifyEqual(testCase, numel(labels), 2);
verifyEqual(testCase, labels(1).kind, 'produced');
verifyEqual(testCase, labels(1).label, "");
clear cleanup
end

function testTxtMissingKindThrows(testCase)
% a txt file without labels and no default should raise an error.
file = writeTempTxt({sprintf('0.00\t0.10')});
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelKind');
clear cleanup
end

function testUnsupportedExtension(testCase)
% create a bogus csv file to check the extension guard.
file = [tempname, '.csv'];
fid = fopen(file, 'w');
assert(fid ~= -1);
fprintf(fid, 'dummy');
fclose(fid);
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_labels(file), 'glm:UnsupportedLabelFormat');
clear cleanup
end

function testMissingField(testCase)
% drop the offset field to trigger schema validation.
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced'}, ...
        'onset', {0.1}, ...
        'label', {'x'})));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelsStruct');
clear cleanup
end

function testInvalidKind(testCase)
% provide an unsupported kind label to exercise the guard.
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'other'}, ...
        'onset', {0.0}, ...
        'offset', {0.1}, ...
        'label', {'unknown'})));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelKind');
clear cleanup
end

function testInvalidTimes(testCase)
% capture a negative duration event to test time validation.
file = writeTempMat(struct( ...
    'events', struct( ...
        'kind', {'produced'}, ...
        'onset', {0.2}, ...
        'offset', {0.1}, ...
        'label', {'bad'})));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_labels(file), 'glm:InvalidLabelTimes');
clear cleanup
end

function file = writeTempMat(data)
% helper to stash an events struct as a temporary mat file.
file = [tempname, '.mat'];
save(file, '-struct', 'data');
end

function file = writeTempTxt(lines)
% helper to write an audacity-style txt file for testing.
if isstring(lines)
    lines = arrayfun(@char, lines(:), 'UniformOutput', false);
elseif ischar(lines)
    lines = {lines};
elseif iscell(lines)
    % assume caller provided cell array of char rows
else
    error('Invalid test input for writeTempTxt');
end
file = [tempname, '.txt'];
fid = fopen(file, 'w');
assert(fid ~= -1);
for ii = 1:numel(lines)
    fprintf(fid, '%s\n', lines{ii});
end
fclose(fid);
end

function deleteIfExists(file)
% clean up any temporary files emitted during testing.
if exist(file, 'file') == 2
    delete(file);
end
end
