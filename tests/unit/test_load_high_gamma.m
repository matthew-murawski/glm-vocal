function tests = test_load_high_gamma
%TEST_LOAD_HIGH_GAMMA Unit tests for load_high_gamma.
%   tests = TEST_LOAD_HIGH_GAMMA() returns function-based tests.

tests = functiontests(localfunctions);
end

function testLoadValidSingleChannel(testCase) %#ok<INUSD>
% test successful load of valid single-channel file
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 1));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file);

assert(size(hg.power, 1) == 100);
assert(size(hg.power, 2) == 1);
assert(hg.fs == 50);
assert(all(isfinite(hg.power)));
clear cleanup
end

function testSingleChannelSelection(testCase) %#ok<INUSD>
% test selecting a specific channel by ID
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3, ...
    'channel_ids', [2, 5, 7]));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file, 'channel_id', 5);

assert(size(hg.power, 1) == 100);
assert(hg.channel_id == 5);
clear cleanup
end

function testMultiChannelReductionMean(testCase) %#ok<INUSD>
% test multi-channel reduction using mean
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file, 'reduce', 'mean');

assert(size(hg.power, 1) == 100);
assert(strcmp(hg.channel_id, 'reduced'));
assert(all(isfinite(hg.power)));
clear cleanup
end

function testMultiChannelReductionMedian(testCase) %#ok<INUSD>
% test multi-channel reduction using median
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file, 'reduce', 'median');

assert(size(hg.power, 1) == 100);
assert(strcmp(hg.channel_id, 'reduced'));
clear cleanup
end

function testTimeVectorFromMetadata(testCase) %#ok<INUSD>
% test time vector constructed from metadata.timestamps
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.5, ...
    'baseline', 10.0, ...
    'nChannels', 1, ...
    'include_metadata_timestamps', true));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file);

expected_t = 0.5 + (0:99)' / 50;
assert(max(abs(hg.t - expected_t)) < 1e-10);
assert(all(diff(hg.t) > 0));
clear cleanup
end

function testTimeVectorFromT0Fs(testCase) %#ok<INUSD>
% test time vector constructed from lfp_t0 and lfp_fs when metadata.timestamps absent
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 1.0, ...
    'baseline', 10.0, ...
    'nChannels', 1, ...
    'include_metadata_timestamps', false));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file);

expected_t = 1.0 + (0:99)' / 50;
assert(max(abs(hg.t - expected_t)) < 1e-10);
clear cleanup
end

function testMissingFieldError(testCase)
% test error when required field is missing
file = [tempname '.mat'];
bad_data = struct();
save(file, 'bad_data');
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file), 'glm:InvalidHGStruct');
clear cleanup
end

function testMissingSubfieldError(testCase)
% test error when required subfield is missing
file = [tempname '.mat'];
lfp_data = struct();
lfp_data.lfp_data = ones(100, 1);
lfp_data.lfp_t0 = 0;
lfp_data.channel_ids = 1;
% missing lfp_fs
save(file, 'lfp_data');
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file), 'glm:InvalidHGStruct');
clear cleanup
end

function testInvalidDataTypeError(testCase)
% test error when data has wrong type
file = [tempname '.mat'];
lfp_data = struct();
lfp_data.lfp_data = single(ones(100, 1));  % wrong type
lfp_data.lfp_fs = 100;
lfp_data.lfp_t0 = 0;
lfp_data.channel_ids = 1;
save(file, 'lfp_data');
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file), 'glm:InvalidHGData');
clear cleanup
end

function testNegativeValuesWarning(testCase)
% test warning when negative values are present
file = [tempname '.mat'];
lfp_data = struct();
lfp_data.lfp_data = [10; 8; -2; 5; -1; 12];
lfp_data.lfp_fs = 100;
lfp_data.lfp_t0 = 0;
lfp_data.channel_ids = 1;
save(file, 'lfp_data');
cleanup = onCleanup(@() deleteIfExists(file));

verifyWarning(testCase, @() load_high_gamma(file), 'glm:NegativeValues');

% verify they are set to zero
hg = load_high_gamma(file);
assert(all(hg.power >= 0));
clear cleanup
end

function testAllZeroDataValidation(testCase) %#ok<INUSD>
% test that validation catches all-zero data
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 0.0, ...
    'nChannels', 1));
cleanup = onCleanup(@() deleteIfExists(file));

hg = load_high_gamma(file);
validation = validate_inputs_hg(hg);

assert(~validation.valid);
assert(~isempty(validation.errors));
clear cleanup
end

function testInvalidChannelIdError(testCase)
% test error when specified channel_id not found
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3, ...
    'channel_ids', [2, 5, 7]));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file, 'channel_id', 999), 'glm:InvalidInput');
clear cleanup
end

function testConflictingParametersError(testCase)
% test error when both channel_id and reduce are specified
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3));
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file, 'channel_id', 2, 'reduce', 'mean'), ...
    'glm:InvalidInput');
clear cleanup
end

function testFileNotFoundError(testCase)
% test error when file does not exist
verifyError(testCase, @() load_high_gamma('nonexistent_file.mat'), 'glm:FileNotFound');
end

function testDefaultChannelWarning(testCase)
% test warning when multiple channels but no selection
file = writeTempHGMat(struct(...
    'nSamples', 100, ...
    'fs', 50, ...
    't0', 0.0, ...
    'baseline', 10.0, ...
    'nChannels', 3));
cleanup = onCleanup(@() deleteIfExists(file));

verifyWarning(testCase, @() load_high_gamma(file), 'glm:MultipleChannels');
clear cleanup
end

function testInvalidFsError(testCase)
% test error when fs is invalid
file = [tempname '.mat'];
lfp_data = struct();
lfp_data.lfp_data = ones(100, 1);
lfp_data.lfp_fs = -10;  % negative fs
lfp_data.lfp_t0 = 0;
lfp_data.channel_ids = 1;
save(file, 'lfp_data');
cleanup = onCleanup(@() deleteIfExists(file));

verifyError(testCase, @() load_high_gamma(file), 'glm:InvalidHGData');
clear cleanup
end

%% helper functions

function file = writeTempHGMat(opts)
% helper to create a temporary high gamma MAT file
file = [tempname '.mat'];

nSamples = opts.nSamples;
fs = opts.fs;
t0 = opts.t0;
baseline = opts.baseline;

if isfield(opts, 'nChannels')
    nChannels = opts.nChannels;
else
    nChannels = 1;
end

if isfield(opts, 'channel_ids')
    channel_ids = opts.channel_ids(:);
else
    channel_ids = (1:nChannels)';
end

if isfield(opts, 'include_metadata_timestamps')
    include_timestamps = opts.include_metadata_timestamps;
else
    include_timestamps = true;
end

% create power matrix
power_matrix = baseline * ones(nSamples, nChannels);

% build lfp_data struct
lfp_data = struct();
lfp_data.lfp_data = power_matrix;
lfp_data.lfp_fs = fs;
lfp_data.lfp_t0 = t0;
lfp_data.channel_ids = channel_ids;

% build metadata struct
if include_timestamps
    metadata = struct();
    metadata.timestamps = t0 + (0:nSamples-1)' / fs;
    metadata.session_id = 'test_session';
    save(file, 'lfp_data', 'metadata');
else
    save(file, 'lfp_data');
end
end

function deleteIfExists(file)
% helper to safely delete a file if it exists
if exist(file, 'file')
    delete(file);
end
end
