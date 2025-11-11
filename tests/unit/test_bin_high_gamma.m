function tests = test_bin_high_gamma
%TEST_BIN_HIGH_GAMMA Unit tests for bin_high_gamma.
%   tests = TEST_BIN_HIGH_GAMMA() returns function-based tests.

tests = functiontests(localfunctions);
end

function testExactBinMatch(testCase) %#ok<INUSD>
% test when fs = 1/dt (no downsampling or upsampling)
dt = 0.01;
fs = 100;  % 1/dt
duration = 1.0;

% create high gamma data
nSamples = round(duration * fs);
hg = struct();
hg.power = 10 * ones(nSamples, 1);
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

% create stim with matching resolution
nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

% create config
cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

% bin the data
binned = bin_high_gamma(hg, stim, cfg);

% verify output
assert(length(binned.power) == length(stim.t));
assert(all(abs(binned.power - 10) < 1e-10));
assert(binned.n_empty_bins == 0);
end

function testDownsamplingMean(testCase) %#ok<INUSD>
% test downsampling with mean aggregation
dt = 0.1;  % 10 Hz bins
fs = 100;  % 100 Hz sampling
duration = 1.0;

% create high gamma data with known pattern
nSamples = round(duration * fs);
hg = struct();
hg.power = mod((0:nSamples-1)', 10) + 1;  % values 1-10 repeating
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

% create stim
nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

% create config
cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

% bin the data
binned = bin_high_gamma(hg, stim, cfg);

% verify output length
assert(length(binned.power) == length(stim.t));
assert(strcmp(binned.method, 'mean'));
end

function testDownsamplingMedian(testCase) %#ok<INUSD>
% test downsampling with median aggregation
dt = 0.1;
fs = 100;
duration = 1.0;

nSamples = round(duration * fs);
hg = struct();
hg.power = (1:nSamples)';
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'median';

binned = bin_high_gamma(hg, stim, cfg);

assert(length(binned.power) == length(stim.t));
assert(strcmp(binned.method, 'median'));
end

function testDownsamplingRMS(testCase) %#ok<INUSD>
% test downsampling with rms aggregation
dt = 0.1;
fs = 100;
duration = 1.0;

nSamples = round(duration * fs);
hg = struct();
hg.power = 5 * ones(nSamples, 1);
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'rms';

binned = bin_high_gamma(hg, stim, cfg);

% for constant values, rms should equal the value
assert(all(abs(binned.power - 5) < 1e-10));
end

function testDownsamplingMax(testCase) %#ok<INUSD>
% test downsampling with max aggregation
dt = 0.1;
fs = 100;
duration = 1.0;

nSamples = round(duration * fs);
hg = struct();
hg.power = rand(nSamples, 1) * 5 + 1;  % random between 1 and 6
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'max';

binned = bin_high_gamma(hg, stim, cfg);

assert(length(binned.power) == length(stim.t));
assert(strcmp(binned.method, 'max'));
% max should be >= original values
assert(all(binned.power >= 1) && all(binned.power <= 6));
end

function testUpsamplingWarning(testCase)
% test that upsampling triggers warning
dt = 0.001;  % 1 ms bins = 1000 Hz
fs = 50;     % only 50 Hz sampling
duration = 0.5;

nSamples = round(duration * fs);
hg = struct();
hg.power = 10 * ones(nSamples, 1);
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

binned = bin_high_gamma(hg, stim, cfg);

% should have warning about upsampling
assert(~isempty(binned.warnings));
assert(any(contains(binned.warnings, 'lower than bin rate')));
end

function testEmptyBinInterpolation(testCase) %#ok<INUSD>
% test interpolation of empty bins
dt = 0.01;
duration = 1.0;

% create sparse high gamma data with gaps
nSamples = 20;
hg = struct();
hg.power = [10; 12; 14; 16; 18; 20; 22; 24; 26; 28; 30; 32; 34; 36; 38; 40; 42; 44; 46; 48];
hg.fs = 10;  % low sampling rate
hg.t = (0:nSamples-1)' * 0.1;  % samples at 0, 0.1, 0.2, ...

% create fine-grained stim
nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

binned = bin_high_gamma(hg, stim, cfg);

% should have empty bins that were interpolated
assert(binned.n_empty_bins > 0);
assert(~isempty(binned.warnings));
% interpolated values should be between min and max
assert(all(binned.power >= min(hg.power)));
assert(all(binned.power <= max(hg.power)));
end

function testSingleBin(testCase) %#ok<INUSD>
% test edge case with single bin
dt = 1.0;
duration = 1.0;
fs = 100;

nSamples = round(duration * fs);
hg = struct();
hg.power = 15 * ones(nSamples, 1);
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

stim = struct();
stim.t = 0;  % single bin

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

binned = bin_high_gamma(hg, stim, cfg);

assert(length(binned.power) == 1);
assert(abs(binned.power - 15) < 1e-10);
end

function testOutputLength(testCase) %#ok<INUSD>
% test that output always matches stim.t length
dt = 0.01;
duration = 2.5;
fs = 200;

nSamples = round(duration * fs);
hg = struct();
hg.power = rand(nSamples, 1) * 10;
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

binned = bin_high_gamma(hg, stim, cfg);

assert(length(binned.power) == length(stim.t));
assert(length(binned.power) == nBins);
end

function testNonNegativity(testCase) %#ok<INUSD>
% test that output is always non-negative
dt = 0.01;
duration = 1.0;
fs = 100;

nSamples = round(duration * fs);
hg = struct();
hg.power = abs(randn(nSamples, 1)) * 5 + 2;  % all positive
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

binned = bin_high_gamma(hg, stim, cfg);

assert(all(binned.power >= 0));
end

function testMethodsDiffer(testCase) %#ok<INUSD>
% test that different methods produce different results (for non-constant data)
dt = 0.05;
duration = 1.0;
fs = 100;

nSamples = round(duration * fs);
hg = struct();
% create data with variability
rng(42);  % reproducible random
hg.power = abs(randn(nSamples, 1)) * 3 + 5;
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

% bin with different methods
methods = {'mean', 'median', 'rms', 'max'};
results = cell(length(methods), 1);

for i = 1:length(methods)
    cfg = struct();
    cfg.preprocess = struct();
    cfg.preprocess.hg_binning_method = methods{i};
    results{i} = bin_high_gamma(hg, stim, cfg);
end

% results should differ for variable data
% mean vs median
assert(~isequal(results{1}.power, results{2}.power));
% mean vs rms
assert(~isequal(results{1}.power, results{3}.power));
% max should be >= mean
assert(all(results{4}.power >= results{1}.power - 1e-10));
end

function testInvalidMethod(testCase)
% test error for invalid aggregation method
dt = 0.01;
fs = 100;
duration = 1.0;

nSamples = round(duration * fs);
hg = struct();
hg.power = ones(nSamples, 1) * 10;
hg.fs = fs;
hg.t = (0:nSamples-1)' / fs;

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'invalid_method';

verifyError(testCase, @() bin_high_gamma(hg, stim, cfg), 'glm:InvalidInput');
end

function testMissingInputFields(testCase)
% test error when required fields are missing
dt = 0.01;
fs = 100;
duration = 1.0;

nSamples = round(duration * fs);
hg = struct();
hg.power = ones(nSamples, 1) * 10;
% missing fs and t

nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;

cfg = struct();
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

verifyError(testCase, @() bin_high_gamma(hg, stim, cfg), 'glm:InvalidInput');
end
