%% demo_step6_first_integration
% first complete end-to-end integration test - checkpoint 1

clr; clc;

fprintf('\n');
fprintf('========================================================================\n');
fprintf('                  FIRST END-TO-END INTEGRATION TEST                    \n');
fprintf('                         CHECKPOINT 1                                  \n');
fprintf('========================================================================\n');
fprintf('\n');

%% section 1: load synthetic session
fprintf('STEP 1: Loading synthetic single-event session...\n');
fprintf('------------------------------------------------------------------------\n');

data = generate_single_event_hg();

fprintf('Session parameters:\n');
fprintf('  Duration: %.1f s\n', data.ground_truth.duration_s);
fprintf('  Bin width: %.4f s (%.0f Hz)\n', data.cfg.dt, 1/data.cfg.dt);
fprintf('  Time bins: %d\n', length(data.stim.t));
fprintf('  Baseline power: %.1f\n', data.ground_truth.baseline);
fprintf('  Noise std: %.2f\n', data.ground_truth.noise_std);
fprintf('\n');

fprintf('Event:\n');
fprintf('  Type: %s\n', data.events.kind);
fprintf('  Time: %.2f s\n', data.ground_truth.event_time);
fprintf('  Duration: %.2f s\n', data.events.t_off - data.events.t_on);
fprintf('\n');

fprintf('True kernel:\n');
fprintf('  Peak amplitude: %.2f\n', data.ground_truth.kernel_peak);
fprintf('  Peak delay: %.3f s\n', data.ground_truth.kernel_delay);
fprintf('  Width (sigma): %.3f s\n', data.ground_truth.kernel_width);
fprintf('\n');

%% section 2: bin high gamma
fprintf('STEP 2: Binning high gamma to GLM time base...\n');
fprintf('------------------------------------------------------------------------\n');

hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

fprintf('Binning results:\n');
fprintf('  Method: %s\n', hg_binned.method);
fprintf('  Output bins: %d\n', length(hg_binned.power));
fprintf('  Power range: [%.2f, %.2f]\n', min(hg_binned.power), max(hg_binned.power));
fprintf('  Empty bins: %d\n', hg_binned.n_empty_bins);
if ~isempty(hg_binned.warnings)
    fprintf('  Warnings: %d\n', length(hg_binned.warnings));
end
fprintf('\n');

%% section 3: assemble design matrix
fprintf('STEP 3: Assembling design matrix...\n');
fprintf('------------------------------------------------------------------------\n');

[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fprintf('Design matrix:\n');
fprintf('  Dimensions: [%d × %d]\n', size(X, 1), size(X, 2));
fprintf('  Sparsity: %.1f%%\n', 100 * (1 - nnz(X) / numel(X)));
fprintf('  Response range: [%.2f, %.2f]\n', min(y), max(y));
fprintf('\n');

fprintf('Regressors:\n');
if isfield(colmap, 'intercept')
    fprintf('  - intercept: 1 column\n');
end
if isfield(colmap, 'heard_fields')
    for i = 1:length(colmap.heard_fields)
        field = colmap.heard_fields{i};
        fprintf('  - %s: %d columns\n', field, length(colmap.(field).cols));
    end
end
if isfield(colmap, 'produced_fields')
    for i = 1:length(colmap.produced_fields)
        field = colmap.produced_fields{i};
        fprintf('  - %s: %d columns\n', field, length(colmap.(field).cols));
    end
end
if isfield(colmap, 'states')
    fprintf('  - states: %d columns\n', length(colmap.states.cols));
end
fprintf('\n');

%% section 4: fit model
fprintf('STEP 4: Fitting GLM...\n');
fprintf('------------------------------------------------------------------------\n');

fit = fit_glm_map_hg(X, y, data.cfg);

fprintf('Optimization:\n');
fprintf('  Converged: %s\n', mat2str(fit.converged));
fprintf('  Iterations: %d\n', fit.output.iterations);
fprintf('  Final NLL: %.6f\n', fit.nll);
fprintf('  Gradient norm: %.6e\n', fit.grad_norm);
fprintf('\n');

% compute metrics
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
corr_val = corr(y, fit.mu);

fprintf('Fit quality:\n');
fprintf('  R²: %.4f\n', r_squared);
fprintf('  Correlation: %.4f\n', corr_val);
fprintf('  RMSE: %.4f\n', sqrt(mean((y - fit.mu).^2)));
fprintf('\n');

%% section 5: extract kernels
fprintf('STEP 5: Extracting kernels...\n');
fprintf('------------------------------------------------------------------------\n');

kernels = extract_kernels(fit.w, colmap, data.stim);

fprintf('Extracted regressors:\n');
regressor_types = fieldnames(kernels);
for i = 1:length(regressor_types)
    field = regressor_types{i};
    if isstruct(kernels.(field)) && isfield(kernels.(field), 'weights')
        fprintf('  - %s: %d weights\n', field, length(kernels.(field).weights));
    else
        fprintf('  - %s: scalar\n', field);
    end
end
fprintf('\n');

%% section 6: kernel recovery analysis
fprintf('KERNEL RECOVERY ANALYSIS:\n');
fprintf('------------------------------------------------------------------------\n');

fitted_kernel = kernels.heard_any.weights;
fitted_times = kernels.heard_any.times;
true_kernel = data.ground_truth.kernel_response;
true_times = data.ground_truth.kernel_times;

% correlation
kernel_corr = corr(fitted_kernel, true_kernel);
fprintf('  Correlation (fitted vs true): %.4f\n', kernel_corr);

% peak location
[fitted_peak, fitted_peak_idx] = max(fitted_kernel);
[true_peak, true_peak_idx] = max(true_kernel);
fitted_peak_time = fitted_times(fitted_peak_idx);
true_peak_time = true_times(true_peak_idx);
peak_time_error = abs(fitted_peak_time - true_peak_time);

fprintf('  Peak location:\n');
fprintf('    True: %.3f s\n', true_peak_time);
fprintf('    Fitted: %.3f s\n', fitted_peak_time);
fprintf('    Error: %.4f s (%.1f bins)\n', peak_time_error, peak_time_error/data.cfg.dt);

% peak magnitude
peak_mag_ratio = fitted_peak / true_peak;
peak_mag_error = abs(fitted_peak - true_peak) / true_peak;

fprintf('  Peak magnitude:\n');
fprintf('    True: %.3f\n', true_peak);
fprintf('    Fitted: %.3f\n', fitted_peak);
fprintf('    Ratio: %.3f\n', peak_mag_ratio);
fprintf('    Relative error: %.1f%%\n', peak_mag_error * 100);

fprintf('\n');

%% section 7: comprehensive visualization
fprintf('CREATING VISUALIZATION...\n');
fprintf('------------------------------------------------------------------------\n');

figure('Position', [100, 100, 1400, 900]);

% subplot 1: observed vs predicted power trace
subplot(2, 2, 1);
hold on;

% plot traces
plot(data.stim.t, y, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Observed');
plot(data.stim.t, fit.mu, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Predicted');

% mark event
event_time = data.ground_truth.event_time;
yl = ylim;
plot([event_time, event_time], yl, 'b--', 'LineWidth', 2, 'DisplayName', 'Event');

title(sprintf('Observed vs Predicted Power Trace (R^2 = %.3f)', r_squared));
xlabel('Time (s)');
ylabel('Power');
legend('Location', 'best');
grid on;
hold off;

% subplot 2: scatter plot
subplot(2, 2, 2);
scatter(y, fit.mu, 20, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;

% diagonal line
lims = [min([y; fit.mu]), max([y; fit.mu])];
plot(lims, lims, 'r--', 'LineWidth', 2);

title(sprintf('Observed vs Predicted (r = %.3f)', corr_val));
xlabel('Observed Power');
ylabel('Predicted Power');
axis equal;
grid on;
hold off;

% subplot 3: fitted vs true kernel
subplot(2, 2, 3);
hold on;

plot(true_times, true_kernel, 'b-', 'LineWidth', 2.5, 'DisplayName', 'True kernel');
plot(fitted_times, fitted_kernel, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted kernel');

% mark peaks
plot(true_peak_time, true_peak, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(fitted_peak_time, fitted_peak, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

title(sprintf('Kernel Recovery (r = %.3f)', kernel_corr));
xlabel('Time from event (s)');
ylabel('Response amplitude');
legend('Location', 'best');
grid on;
hold off;

% subplot 4: residuals over time
subplot(2, 2, 4);
residuals = y - fit.mu;
plot(data.stim.t, residuals, 'k-', 'LineWidth', 0.5);
hold on;

% zero line
plot(data.stim.t([1 end]), [0 0], 'r--', 'LineWidth', 1.5);

% mark event
yl = ylim;
plot([event_time, event_time], yl, 'b--', 'LineWidth', 1);

title(sprintf('Residuals (mean=%.3f, std=%.3f)', mean(residuals), std(residuals)));
xlabel('Time (s)');
ylabel('Residual');
grid on;
hold off;

fprintf('Visualization created.\n\n');

%% section 8: checkpoint assessment
fprintf('========================================================================\n');
fprintf('                      CHECKPOINT 1 ASSESSMENT                          \n');
fprintf('========================================================================\n');
fprintf('\n');

% evaluate all criteria
criteria_met = true;

fprintf('Convergence:\n');
if fit.converged
    fprintf('  [✓] Converged in %d iterations\n', fit.output.iterations);
else
    fprintf('  [✗] Failed to converge\n');
    criteria_met = false;
end
fprintf('  Gradient norm: %.6e\n', fit.grad_norm);
fprintf('\n');

fprintf('Kernel Recovery:\n');
if kernel_corr > 0.8
    fprintf('  [✓] Correlation: %.4f (> 0.8)\n', kernel_corr);
else
    fprintf('  [✗] Correlation: %.4f (<= 0.8)\n', kernel_corr);
    criteria_met = false;
end

if peak_time_error <= data.cfg.dt
    fprintf('  [✓] Peak location error: %.4f s (<= %.4f s)\n', ...
        peak_time_error, data.cfg.dt);
else
    fprintf('  [✗] Peak location error: %.4f s (> %.4f s)\n', ...
        peak_time_error, data.cfg.dt);
    criteria_met = false;
end

if peak_mag_error < 0.3
    fprintf('  [✓] Peak magnitude error: %.1f%% (< 30%%)\n', peak_mag_error * 100);
else
    fprintf('  [✗] Peak magnitude error: %.1f%% (>= 30%%)\n', peak_mag_error * 100);
    criteria_met = false;
end
fprintf('\n');

fprintf('Prediction Quality:\n');
if r_squared > 0.5
    fprintf('  [✓] R²: %.4f (> 0.5)\n', r_squared);
else
    fprintf('  [✗] R²: %.4f (<= 0.5)\n', r_squared);
    criteria_met = false;
end
fprintf('\n');

% final verdict
fprintf('========================================================================\n');
if criteria_met
    fprintf('                                                                        \n');
    fprintf('    ✓ ✓ ✓ CHECKPOINT 1: First end-to-end integration successful! ✓ ✓ ✓ \n');
    fprintf('                                                                        \n');
    fprintf('    All components working together:\n');
    fprintf('      • Data loading and binning\n');
    fprintf('      • Design matrix assembly\n');
    fprintf('      • GLM fitting with convergence\n');
    fprintf('      • Kernel extraction and recovery\n');
    fprintf('                                                                        \n');
    fprintf('    The pipeline is validated and ready for more complex scenarios.\n');
    fprintf('                                                                        \n');
else
    fprintf('                                                                        \n');
    fprintf('    [!] Some criteria not met. Review diagnostics above.\n');
    fprintf('                                                                        \n');
end
fprintf('========================================================================\n');
fprintf('\n');
