function tests = test_io_pipeline
%TEST_IO_PIPELINE End-to-end test for I/O pipeline (Step 1)

tests = functiontests(localfunctions);
end

function testFullIOPipeline(testCase) %#ok<INUSD>
% full end-to-end test: generate -> load -> validate

% add paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'tests', 'synthetic')));

% create temporary file
tempDir = tempname;
mkdir(tempDir);
testFilepath = fullfile(tempDir, 'test_session.mat');
cleanup = onCleanup(@() rmdir(tempDir, 's'));

% step 1: generate synthetic data
filepath = generate_simple_hg_session(...
    'duration', 5, ...
    'fs', 50, ...
    'baseline', 8.0, ...
    'n_channels', 1, ...
    'filepath', testFilepath);

assert(strcmp(filepath, testFilepath));
assert(exist(filepath, 'file') == 2);

% step 2: load the data
hg = load_high_gamma(filepath);

% verify output structure
assert(isstruct(hg));
requiredFields = {'power', 'fs', 't', 'session_id', 'channel_id'};
for idx = 1:numel(requiredFields)
    assert(isfield(hg, requiredFields{idx}), ...
        sprintf('Missing field: %s', requiredFields{idx}));
end

% verify dimensions and types
nSamples = round(5 * 50);  % duration * fs
assert(isequal(size(hg.power), [nSamples, 1]));
assert(isa(hg.power, 'double'));
assert(hg.fs == 50);
assert(isequal(size(hg.t), [nSamples, 1]));

% verify data properties
assert(all(isfinite(hg.power)));
assert(all(hg.power >= 0));
assert(all(diff(hg.t) > 0));

% verify baseline value (should be constant at 8.0)
assert(abs(mean(hg.power) - 8.0) < 1e-10);
assert(std(hg.power) < 1e-10);

% step 3: validate the data
validation = validate_inputs_hg(hg);

% flat trace will trigger validation warning about no variability
% but it still loads successfully
assert(~isempty(validation.errors) || ~isempty(validation.warnings));

clear cleanup
end

function testMultiChannelPipeline(testCase) %#ok<INUSD>
% test multi-channel generation and reduction

% add paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'tests', 'synthetic')));

% create temporary file
tempDir = tempname;
mkdir(tempDir);
testFilepath = fullfile(tempDir, 'test_multi.mat');
cleanup = onCleanup(@() rmdir(tempDir, 's'));

% generate multi-channel data
filepath = generate_simple_hg_session(...
    'duration', 3, ...
    'fs', 100, ...
    'baseline', 12.0, ...
    'n_channels', 3, ...
    'filepath', testFilepath);

% test mean reduction
hg_mean = load_high_gamma(filepath, 'reduce', 'mean');
assert(strcmp(hg_mean.channel_id, 'reduced'));
assert(abs(mean(hg_mean.power) - 12.0) < 1e-10);

% test median reduction
hg_median = load_high_gamma(filepath, 'reduce', 'median');
assert(strcmp(hg_median.channel_id, 'reduced'));
assert(abs(median(hg_median.power) - 12.0) < 1e-10);

% verify both methods produce same result for flat data
assert(isequal(hg_mean.power, hg_median.power));

clear cleanup
end

function testTimeVectorConstruction(testCase) %#ok<INUSD>
% test both branches of time vector construction

% add paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'tests', 'synthetic')));

% create temporary file
tempDir = tempname;
mkdir(tempDir);
testFilepath = fullfile(tempDir, 'test_time.mat');
cleanup = onCleanup(@() rmdir(tempDir, 's'));

% generate file (will have metadata.timestamps)
filepath = generate_simple_hg_session(...
    'duration', 2, ...
    'fs', 50, ...
    'baseline', 5.0, ...
    'filepath', testFilepath);

hg = load_high_gamma(filepath);

% verify time vector properties
assert(abs(hg.t(1) - 0) < 1e-10);
expected_dt = 1 / 50;
actual_dt = mean(diff(hg.t));
assert(abs(actual_dt - expected_dt) < 1e-10);

clear cleanup
end

function testValidationCatchesProblems(testCase) %#ok<INUSD>
% test that validation catches common issues

% add paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

% create data with all zeros (will fail validation)
nSamples = 100;
fs = 50;
t = (0:nSamples-1)' / fs;

hg_zeros = struct();
hg_zeros.power = zeros(nSamples, 1);
hg_zeros.fs = fs;
hg_zeros.t = t;

validation_zeros = validate_inputs_hg(hg_zeros);
assert(~validation_zeros.valid);
assert(~isempty(validation_zeros.errors));

% create data with negative values (triggers warning)
hg_negative = struct();
hg_negative.power = [10; 5; -2; 8; 12];
hg_negative.fs = fs;
hg_negative.t = (0:4)' / fs;

validation_neg = validate_inputs_hg(hg_negative);
assert(~isempty(validation_neg.warnings));

% create data with non-monotonic time (fails validation)
hg_bad_time = struct();
hg_bad_time.power = ones(5, 1) * 10;
hg_bad_time.fs = fs;
hg_bad_time.t = [0; 0.02; 0.01; 0.04; 0.05];  % not monotonic

validation_time = validate_inputs_hg(hg_bad_time);
assert(~validation_time.valid);
assert(~isempty(validation_time.errors));
end
