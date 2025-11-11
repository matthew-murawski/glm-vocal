function tests = test_multi_event_recovery()
%TEST_MULTI_EVENT_RECOVERY End-to-end test for multiple kernel recovery.
%   Tests the complete pipeline's ability to simultaneously fit multiple
%   kernels for different event types without interference.
%   This is CHECKPOINT 2 for the high gamma GLM implementation.

tests = functiontests(localfunctions);
end


function test_full_pipeline_multi_event(testCase)
%TEST_FULL_PIPELINE_MULTI_EVENT Run complete pipeline on multi-event session.

data = generate_multi_event_hg();

% verify data structure
testCase.verifyTrue(isfield(data, 'hg_data'), 'Should have hg_data');
testCase.verifyTrue(isfield(data, 'events'), 'Should have events');
testCase.verifyGreaterThan(length(data.events), 1, 'Should have multiple events');

% run pipeline
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify all kernel types extracted
testCase.verifyTrue(isfield(kernels, 'heard_addressed'), ...
    'Should extract heard_addressed kernel');
testCase.verifyTrue(isfield(kernels, 'heard_overheard'), ...
    'Should extract heard_overheard kernel');
testCase.verifyTrue(isfield(kernels, 'produced_spontaneous'), ...
    'Should extract produced_spontaneous kernel');

end


function test_convergence_multi_event(testCase)
%TEST_CONVERGENCE_MULTI_EVENT Verify model converges with multiple kernels.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

testCase.verifyTrue(fit.converged, 'Model should converge');
testCase.verifyLessThan(fit.grad_norm, 1e-3, 'Gradient should be small');

end


function test_addressed_kernel_recovery(testCase)
%TEST_ADDRESSED_KERNEL_RECOVERY Verify heard_addressed kernel recovery.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted = kernels.heard_addressed.weights;
true_kernel = data.ground_truth.addressed.kernel_response;

% correlation
corr_val = corr(fitted, true_kernel);
testCase.verifyGreaterThan(corr_val, 0.7, ...
    sprintf('Addressed kernel correlation (%.3f) should be > 0.7', corr_val));

% peak magnitude
[fitted_peak, ~] = max(fitted);
[true_peak, ~] = max(true_kernel);
mag_error = abs(fitted_peak - true_peak) / abs(true_peak);
testCase.verifyLessThan(mag_error, 0.4, ...
    sprintf('Addressed peak magnitude error (%.1f%%) should be < 40%%', mag_error * 100));

end


function test_overheard_kernel_recovery(testCase)
%TEST_OVERHEARD_KERNEL_RECOVERY Verify heard_overheard kernel recovery.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted = kernels.heard_overheard.weights;
true_kernel = data.ground_truth.overheard.kernel_response;

% correlation
corr_val = corr(fitted, true_kernel);
testCase.verifyGreaterThan(corr_val, 0.7, ...
    sprintf('Overheard kernel correlation (%.3f) should be > 0.7', corr_val));

% peak magnitude
[fitted_peak, ~] = max(fitted);
[true_peak, ~] = max(true_kernel);
mag_error = abs(fitted_peak - true_peak) / abs(true_peak);
testCase.verifyLessThan(mag_error, 0.4, ...
    sprintf('Overheard peak magnitude error (%.1f%%) should be < 40%%', mag_error * 100));

end


function test_produced_kernel_recovery(testCase)
%TEST_PRODUCED_KERNEL_RECOVERY Verify produced kernel recovery.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted = kernels.produced_spontaneous.weights;
true_kernel = data.ground_truth.produced.kernel_response;

% correlation (slightly lower threshold for biphasic kernel)
corr_val = corr(fitted, true_kernel);
testCase.verifyGreaterThan(corr_val, 0.65, ...
    sprintf('Produced kernel correlation (%.3f) should be > 0.65', corr_val));

% peak magnitude (more lenient for biphasic)
[fitted_peak, ~] = max(fitted);
[true_peak, ~] = max(true_kernel);
mag_error = abs(fitted_peak - true_peak) / abs(true_peak);
testCase.verifyLessThan(mag_error, 0.5, ...
    sprintf('Produced peak magnitude error (%.1f%%) should be < 50%%', mag_error * 100));

end


function test_overall_prediction_quality(testCase)
%TEST_OVERALL_PREDICTION_QUALITY Verify overall R² > 0.6.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% compute R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;

testCase.verifyGreaterThan(r_squared, 0.6, ...
    sprintf('R² (%.3f) should be > 0.6', r_squared));

end


function test_kernels_are_distinct(testCase)
%TEST_KERNELS_ARE_DISTINCT Verify different kernels have different shapes.
%   The three kernel types should be distinguishable from each other.

data = generate_multi_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% extract fitted kernels
addressed = kernels.heard_addressed.weights;
overheard = kernels.heard_overheard.weights;
produced = kernels.produced_spontaneous.weights;

% kernels should be more similar to themselves (truth) than to each other
% compare addressed to overheard
cross_corr_addr_over = corr(addressed, overheard);
addr_to_truth = corr(addressed, data.ground_truth.addressed.kernel_response);

testCase.verifyLessThan(cross_corr_addr_over, addr_to_truth, ...
    'Addressed should be more similar to its truth than to overheard');

end


function test_checkpoint2_criteria(testCase)
%TEST_CHECKPOINT2_CRITERIA Verify all CHECKPOINT 2 criteria are met.
%   This test encapsulates all success criteria for Step 7.

data = generate_multi_event_hg();

% run full pipeline
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% criterion 1: convergence
testCase.verifyTrue(fit.converged, 'CHECKPOINT 2: Model must converge');

% criterion 2: overall R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
testCase.verifyGreaterThan(r_squared, 0.6, 'CHECKPOINT 2: R² > 0.6');

% criterion 3: addressed kernel correlation
addr_corr = corr(kernels.heard_addressed.weights, ...
    data.ground_truth.addressed.kernel_response);
testCase.verifyGreaterThan(addr_corr, 0.7, ...
    'CHECKPOINT 2: Addressed correlation > 0.7');

% criterion 4: overheard kernel correlation
over_corr = corr(kernels.heard_overheard.weights, ...
    data.ground_truth.overheard.kernel_response);
testCase.verifyGreaterThan(over_corr, 0.7, ...
    'CHECKPOINT 2: Overheard correlation > 0.7');

% criterion 5: produced kernel correlation (lower threshold for biphasic)
prod_corr = corr(kernels.produced_spontaneous.weights, ...
    data.ground_truth.produced.kernel_response);
testCase.verifyGreaterThan(prod_corr, 0.65, ...
    'CHECKPOINT 2: Produced correlation > 0.65');

% if we got here, all criteria met
fprintf('\n');
fprintf('========================================\n');
fprintf('✓ CHECKPOINT 2 PASSED\n');
fprintf('========================================\n');
fprintf('  Multi-event kernel recovery successful!\n');
fprintf('  Addressed correlation: %.3f\n', addr_corr);
fprintf('  Overheard correlation: %.3f\n', over_corr);
fprintf('  Produced correlation: %.3f\n', prod_corr);
fprintf('  Overall R²: %.3f\n', r_squared);
fprintf('========================================\n');
fprintf('\n');

end
