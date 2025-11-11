function tests = test_extract_kernels()
%TEST_EXTRACT_KERNELS Unit tests for kernel extraction utility.
%   Tests that extract_kernels correctly parses weight vectors using colmap.

tests = functiontests(localfunctions);
end


function test_basic_extraction(testCase)
%TEST_BASIC_EXTRACTION Test with known weight vector and colmap.

% create simple synthetic weight vector and colmap
dt = 0.01;
n_lags = 10;

% build weight vector: intercept + one heard kernel + states
w = [5.0;  % intercept
     randn(n_lags, 1);  % heard_any kernel
     2.0;  % convo state
     -1.0];  % spon state

% build colmap
colmap = struct();
colmap.intercept = struct('cols', 1, 'name', 'intercept');

% heard_any with info
colmap.heard_any = struct();
colmap.heard_any.cols = 2:(n_lags+1);
colmap.heard_any.info = struct();
colmap.heard_any.info.window_s = [0, (n_lags-1)*dt];
colmap.heard_any.info.lag_mode = 'causal';
colmap.heard_any.info.lag_bins = 0:(n_lags-1);
colmap.heard_any.info.lag_times_s = ((0:(n_lags-1)) * dt)';

colmap.heard_fields = {'heard_any'};

% states
colmap.states = struct();
colmap.states.cols = [n_lags+2, n_lags+3];
colmap.states.names = {'convo', 'spon'};

% stim struct
stim = struct('dt', dt);

% extract kernels
kernels = extract_kernels(w, colmap, stim);

% verify intercept
testCase.verifyTrue(isfield(kernels, 'intercept'), 'Should extract intercept');
testCase.verifyEqual(kernels.intercept, 5.0, 'Intercept should match');

% verify heard_any kernel
testCase.verifyTrue(isfield(kernels, 'heard_any'), 'Should extract heard_any');
testCase.verifyEqual(length(kernels.heard_any.weights), n_lags, ...
    'Kernel should have correct length');
testCase.verifyEqual(length(kernels.heard_any.times), n_lags, ...
    'Time vector should have correct length');

% verify time vector spacing
testCase.verifyEqual(kernels.heard_any.times, ((0:(n_lags-1)) * dt)', ...
    'AbsTol', 1e-10, 'Time vector should be correct');

% verify states
testCase.verifyTrue(isfield(kernels, 'states'), 'Should extract states');
testCase.verifyEqual(kernels.states.convo, 2.0, 'Convo state should match');
testCase.verifyEqual(kernels.states.spon, -1.0, 'Spon state should match');

end


function test_multiple_kernel_types(testCase)
%TEST_MULTIPLE_KERNEL_TYPES Test with heard and produced kernels.

dt = 0.01;
n_lags_heard = 5;
n_lags_produced = 8;

% weight vector
w = [10.0;  % intercept
     ones(n_lags_heard, 1);  % heard_addressed
     2*ones(n_lags_produced, 1);  % produced_spontaneous
     1.5;  % convo
     0.5];  % spon

% colmap
colmap = struct();
colmap.intercept = struct('cols', 1);

% heard kernel
colmap.heard_addressed = struct();
colmap.heard_addressed.cols = 2:(n_lags_heard+1);
colmap.heard_addressed.info = struct();
colmap.heard_addressed.info.window_s = [0, (n_lags_heard-1)*dt];
colmap.heard_addressed.info.lag_times_s = ((0:(n_lags_heard-1)) * dt)';
colmap.heard_addressed.info.lag_bins = 0:(n_lags_heard-1);
colmap.heard_fields = {'heard_addressed'};

% produced kernel
colmap.produced_spontaneous = struct();
start_col = n_lags_heard + 2;
colmap.produced_spontaneous.cols = start_col:(start_col + n_lags_produced - 1);
colmap.produced_spontaneous.info = struct();
colmap.produced_spontaneous.info.window_s = [-0.02, (n_lags_produced-2)*dt];
colmap.produced_spontaneous.info.lag_times_s = ((-0.02:dt:(n_lags_produced-2)*dt))';
colmap.produced_spontaneous.info.lag_bins = -2:(n_lags_produced-3);
colmap.produced_fields = {'produced_spontaneous'};

% states
colmap.states = struct();
colmap.states.cols = [start_col + n_lags_produced, start_col + n_lags_produced + 1];
colmap.states.names = {'convo', 'spon'};

stim = struct('dt', dt);

% extract
kernels = extract_kernels(w, colmap, stim);

% verify both kernels present
testCase.verifyTrue(isfield(kernels, 'heard_addressed'), ...
    'Should have heard_addressed');
testCase.verifyTrue(isfield(kernels, 'produced_spontaneous'), ...
    'Should have produced_spontaneous');

% verify lengths
testCase.verifyEqual(length(kernels.heard_addressed.weights), n_lags_heard, ...
    'Heard kernel length correct');
testCase.verifyEqual(length(kernels.produced_spontaneous.weights), n_lags_produced, ...
    'Produced kernel length correct');

% verify values
testCase.verifyEqual(kernels.heard_addressed.weights, ones(n_lags_heard, 1), ...
    'Heard weights correct');
testCase.verifyEqual(kernels.produced_spontaneous.weights, 2*ones(n_lags_produced, 1), ...
    'Produced weights correct');

end


function test_missing_regressors(testCase)
%TEST_MISSING_REGRESSORS Test with empty regressor lists.

dt = 0.01;

% only intercept and states
w = [5.0;  % intercept
     1.0;  % convo
     0.0];  % spon

colmap = struct();
colmap.intercept = struct('cols', 1);
colmap.states = struct('cols', [2, 3], 'names', {{'convo', 'spon'}});

% no heard or produced fields
colmap.heard_fields = {};
colmap.produced_fields = {};

stim = struct('dt', dt);

% extract (should not error)
kernels = extract_kernels(w, colmap, stim);

% verify intercept and states present, but no kernels
testCase.verifyTrue(isfield(kernels, 'intercept'), 'Should have intercept');
testCase.verifyTrue(isfield(kernels, 'states'), 'Should have states');

% verify no heard/produced kernels
testCase.verifyFalse(isfield(kernels, 'heard_any'), 'Should not have heard_any');
testCase.verifyFalse(isfield(kernels, 'produced_spontaneous'), ...
    'Should not have produced_spontaneous');

end


function test_time_vector_correctness(testCase)
%TEST_TIME_VECTOR_CORRECTNESS Verify time vectors have correct spacing.

dt = 0.01;
n_lags = 20;

% causal kernel [0, 0.19s]
w = [1.0; randn(n_lags, 1)];

colmap = struct();
colmap.intercept = struct('cols', 1);
colmap.heard_any = struct();
colmap.heard_any.cols = 2:(n_lags+1);
colmap.heard_any.info = struct();
colmap.heard_any.info.window_s = [0, (n_lags-1)*dt];
colmap.heard_any.info.lag_times_s = ((0:(n_lags-1)) * dt)';
colmap.heard_any.info.lag_bins = 0:(n_lags-1);
colmap.heard_fields = {'heard_any'};

stim = struct('dt', dt);

kernels = extract_kernels(w, colmap, stim);

% verify time vector
times = kernels.heard_any.times;
testCase.verifyEqual(length(times), n_lags, 'Time vector length');
testCase.verifyEqual(times(1), 0, 'First time should be 0');
testCase.verifyEqual(times(end), (n_lags-1)*dt, 'AbsTol', 1e-10, ...
    'Last time should match window');

% verify uniform spacing
diffs = diff(times);
testCase.verifyEqual(diffs, dt*ones(n_lags-1, 1), 'AbsTol', 1e-10, ...
    'Time spacing should be uniform');

end


function test_integration_with_real_fit(testCase)
%TEST_INTEGRATION_WITH_REAL_FIT Test with actual fitted model.

% generate simple synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% extract kernels
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify structure
testCase.verifyTrue(isfield(kernels, 'intercept'), 'Should have intercept');
testCase.verifyTrue(isfield(kernels, 'heard_any'), 'Should have heard_any kernel');

% verify heard kernel has reasonable properties
testCase.verifyGreaterThan(length(kernels.heard_any.weights), 0, ...
    'Heard kernel should have weights');
testCase.verifyEqual(length(kernels.heard_any.weights), ...
    length(kernels.heard_any.times), ...
    'Weights and times should match');

% verify time vector is monotonic
times = kernels.heard_any.times;
testCase.verifyTrue(all(diff(times) > 0), 'Time vector should be increasing');

end


function test_error_handling(testCase)
%TEST_ERROR_HANDLING Test that errors are caught appropriately.

dt = 0.01;
w = randn(10, 1);
colmap = struct();
stim = struct('dt', dt);

% test with non-vector w
testCase.verifyError(@() extract_kernels(randn(3, 3), colmap, stim), ...
    'glm:InvalidInput', 'Should error on non-vector w');

% test with non-struct colmap
testCase.verifyError(@() extract_kernels(w, [], stim), ...
    'glm:InvalidInput', 'Should error on non-struct colmap');

% test with missing dt
testCase.verifyError(@() extract_kernels(w, colmap, struct()), ...
    'glm:InvalidInput', 'Should error on missing dt');

end
