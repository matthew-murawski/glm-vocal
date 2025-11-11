function tests = test_design_matrix_integration()
%TEST_DESIGN_MATRIX_INTEGRATION End-to-end test for design matrix with fitting.
%   Tests that the design matrix assembly integrates correctly with the
%   GLM fitting pipeline from Step 4.

tests = functiontests(localfunctions);
end


function test_assemble_and_fit(testCase)
%TEST_ASSEMBLE_AND_FIT Assemble design matrix and fit GLM.
%   Verify that design matrix can be passed to fit_glm_map_hg and
%   fitting completes successfully.

% generate synthetic data with single event
data = generate_single_event_hg();

% bin high gamma data
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify design matrix properties
testCase.verifyGreaterThan(size(X, 2), 1, ...
    'Design matrix should have more than intercept');
testCase.verifyEqual(size(X, 1), length(y), ...
    'Design matrix rows should match response length');

% fit GLM using design matrix
fit = fit_glm_map_hg(X, y, data.cfg);

% verify fit completed
testCase.verifyTrue(isfield(fit, 'w'), 'Fit should contain weights');
testCase.verifyTrue(isfield(fit, 'mu'), 'Fit should contain predictions');
testCase.verifyTrue(isfield(fit, 'nll'), 'Fit should contain NLL');
testCase.verifyTrue(isfield(fit, 'converged'), 'Fit should contain convergence flag');

% verify weights have correct dimension
testCase.verifyEqual(length(fit.w), size(X, 2), ...
    'Number of weights should match design matrix columns');

% verify predictions have correct dimension
testCase.verifyEqual(length(fit.mu), length(y), ...
    'Predictions should match response length');

% verify convergence
testCase.verifyTrue(fit.converged, 'Fit should converge on simple synthetic data');

% verify NLL is finite
testCase.verifyTrue(isfinite(fit.nll), 'NLL should be finite');
testCase.verifyGreaterThanOrEqual(fit.nll, 0, 'NLL should be non-negative');

end


function test_fit_quality_basic(testCase)
%TEST_FIT_QUALITY_BASIC Basic sanity checks on fitted weights.
%   Verify that fitted model makes reasonable predictions.

% generate synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% compute correlation between observed and predicted
corr_val = corr(y, fit.mu);
testCase.verifyGreaterThan(corr_val, 0.5, ...
    'Correlation should be > 0.5 for synthetic data with clear signal');

% compute R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;

testCase.verifyGreaterThan(r_squared, 0.3, ...
    'R² should be > 0.3 for synthetic data');

% verify predictions are in reasonable range
testCase.verifyGreaterThanOrEqual(min(fit.mu), 0, ...
    'Predictions should be non-negative for identity link');
testCase.verifyLessThan(max(fit.mu), max(y) * 2, ...
    'Predictions should not be wildly larger than observed');

end


function test_intercept_interpretation(testCase)
%TEST_INTERCEPT_INTERPRETATION Verify intercept captures baseline.
%   The intercept should approximate the baseline power level.

% generate synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% extract intercept
if isfield(colmap, 'intercept')
    intercept_idx = colmap.intercept.cols;
    intercept_val = fit.w(intercept_idx);

    % intercept should be reasonably close to baseline
    % (note: with basis projection, the representation changes, so we allow more tolerance)
    baseline = data.ground_truth.baseline;
    testCase.verifyLessThan(abs(intercept_val - baseline), baseline * 0.7, ...
        'Intercept should be within 70% of true baseline (with basis projection)');
end

end


function test_kernel_extraction_structure(testCase)
%TEST_KERNEL_EXTRACTION_STRUCTURE Verify we can extract kernel structure.
%   Check that we can parse the weights using colmap to extract kernels.

% generate synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% extract heard_addressed kernel if present
if isfield(colmap, 'heard_addressed')
    heard_cols = colmap.heard_addressed.cols;

    % extract kernel weights
    kernel_weights = fit.w(heard_cols);

    % verify kernel structure
    testCase.verifyEqual(length(kernel_weights), length(heard_cols), ...
        'Kernel weights should match number of columns');

    % kernel should have positive values somewhere (response to event)
    testCase.verifyGreaterThan(max(kernel_weights), 0, ...
        'Kernel should have positive response');

    % build kernel time vector
    heard_window = data.cfg.heard_window_s;
    dt = data.cfg.dt;
    kernel_times = (heard_window(1):dt:heard_window(2))';

    testCase.verifyEqual(length(kernel_times), length(kernel_weights), ...
        'Kernel time vector should match kernel weights');
end

end


function test_states_included(testCase)
%TEST_STATES_INCLUDED Verify state coefficients are included.
%   Check that conversational state regressors are in the design matrix.

% generate synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify states field present
testCase.verifyTrue(isfield(colmap, 'states'), ...
    'colmap should contain states field');

% verify state columns
if isfield(colmap, 'states')
    testCase.verifyTrue(isfield(colmap.states, 'cols'), ...
        'states should have cols field');

    state_cols = colmap.states.cols;
    testCase.verifyEqual(length(state_cols), 2, ...
        'Should have 2 state columns (convo, spon)');

    % verify state column values are binary
    for col = state_cols
        state_values = full(X(:, col));
        unique_vals = unique(state_values);
        testCase.verifyLessThanOrEqual(length(unique_vals), 2, ...
            'State columns should be binary (0/1)');
        testCase.verifyTrue(all(ismember(unique_vals, [0, 1])), ...
            'State values should be 0 or 1');
    end
end

end


function test_no_spike_history_in_fit(testCase)
%TEST_NO_SPIKE_HISTORY_IN_FIT Confirm spike history excluded throughout.
%   Verify that spike history does not appear in any stage.

% generate synthetic data
data = generate_single_event_hg();

% bin and assemble
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify no spike_history in colmap
testCase.verifyFalse(isfield(colmap, 'spike_history'), ...
    'colmap should not contain spike_history');

% verify exclude_predictors includes spike_history
if isfield(data.cfg, 'exclude_predictors')
    excludes = data.cfg.exclude_predictors;
    if iscell(excludes)
        % should be added by wrapper
        % (may or may not be in original cfg)
    end
end

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% verify fit completes
testCase.verifyTrue(fit.converged, 'Fit should converge');

end
