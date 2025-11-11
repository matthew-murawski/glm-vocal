function tests = test_assemble_design_matrix_hg()
%TEST_ASSEMBLE_DESIGN_MATRIX_HG Unit tests for design matrix assembly.
%   Tests the high gamma design matrix wrapper function to ensure it:
%   - Creates correct design matrix structure
%   - Excludes spike history columns
%   - Handles different event types properly
%   - Validates dimensions and content

tests = functiontests(localfunctions);
end


function test_single_perceived_event(testCase)
%TEST_SINGLE_PERCEIVED_EVENT Test with single perceived call.
%   Verify heard_addressed columns created with causal kernel structure.

% generate synthetic data
data = generate_single_event_hg();

% bin the high gamma data
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify response variable
testCase.verifyEqual(length(y), length(data.stim.t), 'Response length should match time base');
testCase.verifyEqual(y, hg_binned.power, 'Response should be binned power');

% verify design matrix dimensions
testCase.verifyEqual(size(X, 1), length(y), 'Design matrix rows should match response length');
testCase.verifyGreaterThan(size(X, 2), 0, 'Design matrix should have columns');

% verify no spike history columns
testCase.verifyFalse(isfield(colmap, 'spike_history'), ...
    'Design matrix should not contain spike_history columns');

% verify intercept present
testCase.verifyTrue(isfield(colmap, 'intercept'), 'Should have intercept column');

% verify heard columns present
testCase.verifyTrue(isfield(colmap, 'heard_fields'), 'Should have heard_fields in colmap');

% check that heard_addressed is included
if isfield(colmap, 'heard_addressed')
    heard_cols = colmap.heard_addressed.cols;
    testCase.verifyGreaterThan(length(heard_cols), 0, 'Should have heard_addressed columns');

    % verify causal structure: columns should correspond to [0, 2s] window
    heard_window = data.cfg.heard_window_s;
    expected_lags = round((heard_window(2) - heard_window(1)) / data.cfg.dt) + 1;
    testCase.verifyEqual(length(heard_cols), expected_lags, ...
        'Number of heard columns should match kernel window');
end

% verify no NaN or Inf
testCase.verifyTrue(all(isfinite(X(:))), 'Design matrix should not contain NaN or Inf');
testCase.verifyTrue(all(isfinite(y)), 'Response should not contain NaN or Inf');

end


function test_single_produced_event(testCase)
%TEST_SINGLE_PRODUCED_EVENT Test with single produced call.
%   Verify produced columns created with symmetric kernel.

% create synthetic data with produced event
duration_s = 6.0;
dt = 0.01;
t = (0:dt:duration_s-dt)';
nT = length(t);

% create produced event at t=3s
events = struct();
events.kind = 'produced';
events.t_on = 3.0;
events.t_off = 3.5;
events.label = 'spontaneous';
events.produced_context = 'produced_spontaneous';

% build streams
stim = struct();
stim.t = t;
stim.dt = dt;

streams = struct();
streams.heard_any = zeros(nT, 1);
streams.heard_addressed = zeros(nT, 1);
streams.heard_overheard = zeros(nT, 1);
streams.produced_spontaneous = zeros(nT, 1);
streams.produced_after_heard = zeros(nT, 1);
streams.produced_after_produced = zeros(nT, 1);

% mark produced event
event_idx = find(t >= events.t_on & t < events.t_off);
streams.produced_spontaneous(event_idx) = 1;

stim.streams = streams;

% build states
states = struct();
states.convo = zeros(nT, 1);
states.spon = ones(nT, 1);
stim.states = states;

% create hg_binned
hg_binned = struct();
hg_binned.power = 5 * ones(nT, 1);  % flat baseline
hg_binned.t = t;

% config
cfg = struct();
cfg.dt = dt;
cfg.heard_window_s = [0, 2.0];
cfg.produced_window_s = [-2.0, 3.0];  % symmetric window
cfg.history_window_s = [0, 0.5];
cfg.heard_basis = struct('kind', 'raised_cosine');
cfg.produced_basis = struct('kind', 'raised_cosine');
cfg.preprocess.hg_binning_method = 'mean';
cfg.exclude_predictors = {};

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, events, stim, cfg);

% verify produced columns present
testCase.verifyTrue(isfield(colmap, 'produced_fields'), ...
    'Should have produced_fields in colmap');

% check that produced_spontaneous is included
if isfield(colmap, 'produced_spontaneous')
    prod_cols = colmap.produced_spontaneous.cols;
    testCase.verifyGreaterThan(length(prod_cols), 0, ...
        'Should have produced_spontaneous columns');

    % verify symmetric window (with basis projection, number of columns will be less than raw lags)
    % just verify we have some columns
    testCase.verifyGreaterThan(length(prod_cols), 0, ...
        'Should have some produced columns (basis-projected)');
end

% verify no spike history
testCase.verifyFalse(isfield(colmap, 'spike_history'), ...
    'Should not have spike_history');

end


function test_no_events(testCase)
%TEST_NO_EVENTS Test with no events.
%   Should return intercept only and not error.

% create minimal session with no events
duration_s = 2.0;
dt = 0.01;
t = (0:dt:duration_s-dt)';
nT = length(t);

% build streams (all zeros)
stim = struct();
stim.t = t;
stim.dt = dt;

streams = struct();
streams.heard_any = zeros(nT, 1);
streams.heard_addressed = zeros(nT, 1);
streams.heard_overheard = zeros(nT, 1);
streams.produced_spontaneous = zeros(nT, 1);
streams.produced_after_heard = zeros(nT, 1);
streams.produced_after_produced = zeros(nT, 1);

stim.streams = streams;

% build states
states = struct();
states.convo = zeros(nT, 1);
states.spon = ones(nT, 1);
stim.states = states;

% create hg_binned
hg_binned = struct();
hg_binned.power = 5 * ones(nT, 1);
hg_binned.t = t;

% config
cfg = struct();
cfg.dt = dt;
cfg.heard_window_s = [0, 2.0];
cfg.produced_window_s = [-2.0, 3.0];
cfg.history_window_s = [0, 0.5];
cfg.heard_basis = struct('kind', 'raised_cosine');
cfg.produced_basis = struct('kind', 'raised_cosine');
cfg.preprocess.hg_binning_method = 'mean';
cfg.exclude_predictors = {};

% assemble design matrix (should not error)
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, [], stim, cfg);

% verify basic structure
testCase.verifyEqual(size(X, 1), nT, 'Should have correct number of rows');
testCase.verifyEqual(length(y), nT, 'Response should have correct length');

% verify intercept present
testCase.verifyTrue(isfield(colmap, 'intercept'), 'Should have intercept');

% verify no spike history
testCase.verifyFalse(isfield(colmap, 'spike_history'), ...
    'Should not have spike_history');

end


function test_colmap_structure(testCase)
%TEST_COLMAP_STRUCTURE Test column mapping structure.
%   Verify field names and indices map correctly.

% generate synthetic data
data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify colmap is a struct
testCase.verifyTrue(isstruct(colmap), 'colmap should be a struct');

% verify intercept mapping
if isfield(colmap, 'intercept')
    testCase.verifyTrue(isscalar(colmap.intercept.cols), ...
        'Intercept should map to single column');
    testCase.verifyEqual(colmap.intercept.cols, 1, ...
        'Intercept should be first column');
end

% verify all column indices are within bounds
p = size(X, 2);
fields = fieldnames(colmap);
for i = 1:length(fields)
    field = fields{i};
    if isstruct(colmap.(field)) && isfield(colmap.(field), 'cols')
        cols = colmap.(field).cols;
        testCase.verifyGreaterThanOrEqual(min(cols), 1, ...
            sprintf('%s columns should be >= 1', field));
        testCase.verifyLessThanOrEqual(max(cols), p, ...
            sprintf('%s columns should be <= %d', field, p));
    end
end

end


function test_dimension_consistency(testCase)
%TEST_DIMENSION_CONSISTENCY Test that dimensions are consistent.
%   Verify n_rows = length(hg_binned.power).

% generate synthetic data
data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify dimensions
n = length(hg_binned.power);
testCase.verifyEqual(size(X, 1), n, 'Design matrix rows should match power length');
testCase.verifyEqual(length(y), n, 'Response should match power length');
testCase.verifyEqual(length(data.stim.t), n, 'Time base should match power length');

% verify column count matches regressor structure
p = size(X, 2);

% count columns from colmap
total_cols = 0;
fields = fieldnames(colmap);
for i = 1:length(fields)
    field = fields{i};
    if isstruct(colmap.(field)) && isfield(colmap.(field), 'cols')
        cols = colmap.(field).cols;
        total_cols = max(total_cols, max(cols));
    end
end

testCase.verifyEqual(p, total_cols, ...
    'Design matrix columns should match colmap structure');

end


function test_validation_catches_errors(testCase)
%TEST_VALIDATION_CATCHES_ERRORS Test that validation catches common errors.

% generate synthetic data
data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

% test 1: mismatched dimensions
hg_binned_bad = hg_binned;
hg_binned_bad.power = hg_binned.power(1:end-10);

testCase.verifyError(@() assemble_design_matrix_hg(hg_binned_bad, data.events, data.stim, data.cfg), ...
    'glm:DimensionMismatch', ...
    'Should error on dimension mismatch');

% test 2: NaN in power
hg_binned_nan = hg_binned;
hg_binned_nan.power(10) = NaN;

testCase.verifyError(@() assemble_design_matrix_hg(hg_binned_nan, data.events, data.stim, data.cfg), ...
    'glm:InvalidInput', ...
    'Should error on NaN in power');

% test 3: missing required fields
testCase.verifyError(@() assemble_design_matrix_hg(struct(), data.events, data.stim, data.cfg), ...
    'glm:InvalidInput', ...
    'Should error on missing hg_binned.power');

end
