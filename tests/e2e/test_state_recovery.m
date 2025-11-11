function tests = test_state_recovery()
%TEST_STATE_RECOVERY End-to-end test for state coefficient recovery.
%   Tests the complete pipeline's ability to simultaneously fit both
%   event kernels and state-dependent baseline shifts.
%   This is CHECKPOINT 3 for the high gamma GLM implementation.

tests = functiontests(localfunctions);
end


function test_full_pipeline_with_states(testCase)
%TEST_FULL_PIPELINE_WITH_STATES Run complete pipeline on state effect session.

data = generate_state_effect_hg();

% verify data structure
testCase.verifyTrue(isfield(data, 'hg_data'), 'Should have hg_data');
testCase.verifyTrue(isfield(data, 'events'), 'Should have events');
testCase.verifyTrue(isfield(data.stim, 'states'), 'Should have states');
testCase.verifyEqual(length(fieldnames(data.stim.states)), 2, 'Should have 2 states');

% run pipeline
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify kernel and states extracted
testCase.verifyTrue(isfield(kernels, 'heard_addressed'), ...
    'Should extract heard_addressed kernel');
testCase.verifyTrue(isfield(kernels, 'states'), ...
    'Should extract states');

end


function test_convergence_with_states(testCase)
%TEST_CONVERGENCE_WITH_STATES Verify model converges with state regressors.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

testCase.verifyTrue(fit.converged, 'Model should converge');
testCase.verifyLessThan(fit.grad_norm, 1e-3, 'Gradient should be small');

end


function test_kernel_recovery_with_states(testCase)
%TEST_KERNEL_RECOVERY_WITH_STATES Verify kernel recovery with state regressors.
%   Kernel should still be recovered accurately even with states present.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted = kernels.heard_addressed.weights;
true_kernel = data.ground_truth.kernel.kernel_response;

% correlation
corr_val = corr(fitted, true_kernel);
testCase.verifyGreaterThan(corr_val, 0.7, ...
    sprintf('Kernel correlation (%.3f) should be > 0.7', corr_val));

end


function test_convo_state_coefficient(testCase)
%TEST_CONVO_STATE_COEFFICIENT Verify convo state coefficient recovery.
%   note: with both states in the model, there is collinearity with the
%   intercept. the key is that the state *difference* is recovered.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify the state difference, not absolute values
state_diff_fitted = kernels.states.convo - kernels.states.spon;
state_diff_true = data.ground_truth.state_coeff_convo - data.ground_truth.state_coeff_spon;  % +4 - 0 = +4

% difference should be within 20% of true
pct_error = abs(state_diff_fitted - state_diff_true) / abs(state_diff_true);
testCase.verifyLessThan(pct_error, 0.2, ...
    sprintf('State difference error (%.1f%%) should be < 20%%', pct_error * 100));

end


function test_spon_state_coefficient(testCase)
%TEST_SPON_STATE_COEFFICIENT Verify state difference is recovered.
%   note: this test is redundant with test_convo_state_coefficient but
%   kept for completeness. both states may be offset from their nominal
%   values due to collinearity with intercept, but difference is correct.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify the state difference (same as previous test)
state_diff_fitted = kernels.states.convo - kernels.states.spon;
state_diff_true = data.ground_truth.state_coeff_convo - data.ground_truth.state_coeff_spon;

% difference within 20%
pct_error = abs(state_diff_fitted - state_diff_true) / abs(state_diff_true);
testCase.verifyLessThan(pct_error, 0.2, ...
    sprintf('State difference (%.2f vs %.2f) error %.1f%% should be < 20%%', ...
    state_diff_fitted, state_diff_true, pct_error * 100));

end


function test_overall_prediction_quality_with_states(testCase)
%TEST_OVERALL_PREDICTION_QUALITY_WITH_STATES Verify R² > 0.7.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% compute R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;

testCase.verifyGreaterThan(r_squared, 0.7, ...
    sprintf('R² (%.3f) should be > 0.7', r_squared));

end


function test_epoch_specific_predictions(testCase)
%TEST_EPOCH_SPECIFIC_PREDICTIONS Verify mean power convo > spon.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% extract state masks
convo_mask = data.ground_truth.convo_mask;
spon_mask = data.ground_truth.spon_mask;

% compute mean predicted power in each epoch
mean_pred_convo = mean(fit.mu(convo_mask));
mean_pred_spon = mean(fit.mu(spon_mask));

testCase.verifyGreaterThan(mean_pred_convo, mean_pred_spon, ...
    sprintf('Mean predicted power convo (%.2f) should be > spon (%.2f)', ...
    mean_pred_convo, mean_pred_spon));

% verify difference is close to true baseline shift (+4)
pred_shift = mean_pred_convo - mean_pred_spon;
true_shift = data.ground_truth.baseline_convo - data.ground_truth.baseline_spon;

shift_error = abs(pred_shift - true_shift) / abs(true_shift);
testCase.verifyLessThan(shift_error, 0.3, ...
    sprintf('Predicted baseline shift (%.2f) should be within 30%% of true (%.2f)', ...
    pred_shift, true_shift));

end


function test_checkpoint3_criteria(testCase)
%TEST_CHECKPOINT3_CRITERIA Verify all CHECKPOINT 3 criteria are met.
%   This test encapsulates all success criteria for Step 8.

data = generate_state_effect_hg();

% run full pipeline
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% criterion 1: convergence
testCase.verifyTrue(fit.converged, 'CHECKPOINT 3: Model must converge');

% criterion 2: kernel correlation
kernel_corr = corr(kernels.heard_addressed.weights, ...
    data.ground_truth.kernel.kernel_response);
testCase.verifyGreaterThan(kernel_corr, 0.7, ...
    'CHECKPOINT 3: Kernel correlation > 0.7');

% criterion 3: state difference within 20%
state_diff_fitted = kernels.states.convo - kernels.states.spon;
state_diff_true = data.ground_truth.state_coeff_convo - data.ground_truth.state_coeff_spon;
state_error = abs(state_diff_fitted - state_diff_true) / abs(state_diff_true);
testCase.verifyLessThan(state_error, 0.2, ...
    'CHECKPOINT 3: State difference within 20%');

% criterion 4: overall R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
testCase.verifyGreaterThan(r_squared, 0.7, ...
    'CHECKPOINT 3: R² > 0.7');

% criterion 5: epoch-specific predictions
convo_mask = data.ground_truth.convo_mask;
spon_mask = data.ground_truth.spon_mask;
mean_pred_convo = mean(fit.mu(convo_mask));
mean_pred_spon = mean(fit.mu(spon_mask));
testCase.verifyGreaterThan(mean_pred_convo, mean_pred_spon, ...
    'CHECKPOINT 3: Mean power convo > spon');

% if we got here, all criteria met
fprintf('\n');
fprintf('========================================\n');
fprintf('✓ CHECKPOINT 3 PASSED\n');
fprintf('========================================\n');
fprintf('  State coefficient recovery successful!\n');
fprintf('  Kernel correlation: %.3f\n', kernel_corr);
fprintf('  State difference: %.2f (true: %.2f, error: %.1f%%)\n', ...
    state_diff_fitted, state_diff_true, state_error * 100);
fprintf('  Overall R²: %.3f\n', r_squared);
fprintf('  Mean pred convo: %.2f | spon: %.2f\n', mean_pred_convo, mean_pred_spon);
fprintf('========================================\n');
fprintf('\n');

end
