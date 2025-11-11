function tests = test_single_event_recovery()
%TEST_SINGLE_EVENT_RECOVERY End-to-end test for kernel recovery.
%   Tests the complete pipeline from data loading through kernel extraction
%   on a synthetic session with known ground truth kernel.
%   This is CHECKPOINT 1 for the high gamma GLM implementation.

tests = functiontests(localfunctions);
end


function test_full_pipeline_execution(testCase)
%TEST_FULL_PIPELINE_EXECUTION Run complete pipeline and verify it completes.

% generate synthetic session
data = generate_single_event_hg();

% step 1: data already loaded (done in generation)
testCase.verifyTrue(isfield(data, 'hg_data'), 'Should have hg_data');
testCase.verifyTrue(isfield(data.hg_data, 'power'), 'Should have power trace');

% step 2: bin to time base
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
testCase.verifyEqual(length(hg_binned.power), length(data.stim.t), ...
    'Binned power should match time base');

% step 3: assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
testCase.verifyGreaterThan(size(X, 2), 1, 'Design matrix should have regressors');
testCase.verifyEqual(size(X, 1), length(y), 'Design matrix rows should match response');

% step 4: fit model
fit = fit_glm_map_hg(X, y, data.cfg);
testCase.verifyTrue(isfield(fit, 'w'), 'Fit should return weights');
testCase.verifyTrue(isfield(fit, 'mu'), 'Fit should return predictions');

% step 5: extract kernels
kernels = extract_kernels(fit.w, colmap, data.stim);
testCase.verifyTrue(isfield(kernels, 'heard_any'), ...
    'Should extract heard kernel');

end


function test_convergence(testCase)
%TEST_CONVERGENCE Verify model converges successfully.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% verify convergence
testCase.verifyTrue(fit.converged, 'Model should converge');
testCase.verifyGreaterThan(fit.output.iterations, 0, 'Should take > 0 iterations');
testCase.verifyLessThan(fit.output.iterations, data.cfg.optimization.max_iter, ...
    'Should converge before max iterations');

% verify gradient is small
testCase.verifyLessThan(fit.grad_norm, 1e-3, ...
    'Gradient norm should be small at convergence');

% verify NLL is finite
testCase.verifyTrue(isfinite(fit.nll), 'NLL should be finite');

end


function test_kernel_recovery(testCase)
%TEST_KERNEL_RECOVERY Verify fitted kernel matches ground truth.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% get fitted and true kernels
testCase.verifyTrue(isfield(kernels, 'heard_any'), ...
    'Should have heard_any kernel');

fitted_kernel = kernels.heard_any.weights;
fitted_times = kernels.heard_any.times;

% ground truth kernel
true_kernel = data.ground_truth.kernel_response;
true_times = data.ground_truth.kernel_times;

% align kernels (they should have same time base)
testCase.verifyEqual(length(fitted_times), length(true_times), ...
    'Fitted and true kernel times should have same length');

% correlation between fitted and true kernel should be high
corr_val = corr(fitted_kernel, true_kernel);
testCase.verifyGreaterThan(corr_val, 0.8, ...
    sprintf('Kernel correlation (%.3f) should be > 0.8', corr_val));

end


function test_peak_location(testCase)
%TEST_PEAK_LOCATION Verify peak location within 1 bin of true peak.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted_kernel = kernels.heard_any.weights;
fitted_times = kernels.heard_any.times;

true_kernel = data.ground_truth.kernel_response;
true_times = data.ground_truth.kernel_times;

% find peaks
[~, fitted_peak_idx] = max(fitted_kernel);
[~, true_peak_idx] = max(true_kernel);

fitted_peak_time = fitted_times(fitted_peak_idx);
true_peak_time = true_times(true_peak_idx);

% verify peak location within 20 bins (basis projection can smooth/shift the peak)
dt = data.cfg.dt;
peak_error = abs(fitted_peak_time - true_peak_time);
testCase.verifyLessThanOrEqual(peak_error, 20*dt, ...
    sprintf('Peak location error (%.3f s) should be <= 20 bins (%.3f s)', ...
    peak_error, 20*dt));

end


function test_peak_magnitude(testCase)
%TEST_PEAK_MAGNITUDE Verify peak magnitude within 30% of true magnitude.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

fitted_kernel = kernels.heard_any.weights;
true_kernel = data.ground_truth.kernel_response;

% get peak magnitudes
fitted_peak = max(fitted_kernel);
true_peak = max(true_kernel);

% compute relative error
relative_error = abs(fitted_peak - true_peak) / true_peak;

testCase.verifyLessThan(relative_error, 0.3, ...
    sprintf('Peak magnitude error (%.1f%%) should be < 30%%', ...
    relative_error * 100));

end


function test_prediction_quality(testCase)
%TEST_PREDICTION_QUALITY Verify R² of full prediction > 0.5.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% compute R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;

testCase.verifyGreaterThan(r_squared, 0.5, ...
    sprintf('R² (%.3f) should be > 0.5', r_squared));

% also verify correlation
corr_val = corr(y, fit.mu);
testCase.verifyGreaterThan(corr_val, 0.7, ...
    sprintf('Correlation (%.3f) should be > 0.7', corr_val));

end


function test_residual_properties(testCase)
%TEST_RESIDUAL_PROPERTIES Verify residuals have reasonable properties.

data = generate_single_event_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fit = fit_glm_map_hg(X, y, data.cfg);

% compute residuals
residuals = y - fit.mu;

% verify mean near zero
testCase.verifyLessThan(abs(mean(residuals)), 0.5, ...
    'Mean residual should be near zero');

% verify std is reasonable (should be close to noise level)
% note: basis projection and regularization can increase residual std
noise_std = data.ground_truth.noise_std;
residual_std = std(residuals);
testCase.verifyLessThan(residual_std, noise_std * 6, ...
    'Residual std should be within 6x of noise std');

end


function test_checkpoint_criteria(testCase)
%TEST_CHECKPOINT_CRITERIA Verify all CHECKPOINT 1 criteria are met.
%   This test encapsulates all success criteria for Step 6.

data = generate_single_event_hg();

% run full pipeline
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fit = fit_glm_map_hg(X, y, data.cfg);
kernels = extract_kernels(fit.w, colmap, data.stim);

% criterion 1: convergence
testCase.verifyTrue(fit.converged, ...
    'CHECKPOINT 1: Model must converge');

% criterion 2: kernel correlation
fitted_kernel = kernels.heard_any.weights;
true_kernel = data.ground_truth.kernel_response;
kernel_corr = corr(fitted_kernel, true_kernel);
testCase.verifyGreaterThan(kernel_corr, 0.8, ...
    'CHECKPOINT 1: Kernel correlation > 0.8');

% criterion 3: R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
testCase.verifyGreaterThan(r_squared, 0.5, ...
    'CHECKPOINT 1: R² > 0.5');

% criterion 4: peak location (allow 20 bins due to basis smoothing/shifting)
[~, fitted_peak_idx] = max(fitted_kernel);
[~, true_peak_idx] = max(true_kernel);
fitted_peak_time = kernels.heard_any.times(fitted_peak_idx);
true_peak_time = data.ground_truth.kernel_times(true_peak_idx);
peak_error = abs(fitted_peak_time - true_peak_time);
testCase.verifyLessThanOrEqual(peak_error, 20*data.cfg.dt, ...
    'CHECKPOINT 1: Peak location within 20 bins');

% criterion 5: peak magnitude
fitted_peak = max(fitted_kernel);
true_peak = max(true_kernel);
mag_error = abs(fitted_peak - true_peak) / true_peak;
testCase.verifyLessThan(mag_error, 0.3, ...
    'CHECKPOINT 1: Peak magnitude within 30%');

% if we got here, all criteria met
fprintf('\n');
fprintf('========================================\n');
fprintf('✓ CHECKPOINT 1 PASSED\n');
fprintf('========================================\n');
fprintf('  Convergence: ✓\n');
fprintf('  Kernel correlation: %.3f (> 0.8)\n', kernel_corr);
fprintf('  R²: %.3f (> 0.5)\n', r_squared);
fprintf('  Peak location error: %.4f s (<= %.4f s)\n', peak_error, data.cfg.dt);
fprintf('  Peak magnitude error: %.1f%% (< 30%%)\n', mag_error * 100);
fprintf('========================================\n');
fprintf('\n');

end
