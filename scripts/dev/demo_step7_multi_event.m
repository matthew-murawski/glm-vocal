%% demo_step7_multi_event
% multi-event kernel recovery test - checkpoint 2

clr; clc;

fprintf('\n');
fprintf('========================================================================\n');
fprintf('                   MULTI-EVENT KERNEL RECOVERY TEST                    \n');
fprintf('                         CHECKPOINT 2                                  \n');
fprintf('========================================================================\n');
fprintf('\n');

%% section 1: load synthetic session
fprintf('LOADING SYNTHETIC MULTI-EVENT SESSION...\n');
fprintf('------------------------------------------------------------------------\n');

data = generate_multi_event_hg();

fprintf('Session parameters:\n');
fprintf('  Duration: %.1f s\n', data.ground_truth.duration_s);
fprintf('  Baseline: %.1f\n', data.ground_truth.baseline);
fprintf('  Noise std: %.2f\n', data.ground_truth.noise_std);
fprintf('\n');

% count events by type
n_addressed = length(data.ground_truth.addressed.event_times);
n_overheard = length(data.ground_truth.overheard.event_times);
n_produced = length(data.ground_truth.produced.event_times);

fprintf('Event summary:\n');
fprintf('  Heard (addressed): %d events at t=[%s]s\n', ...
    n_addressed, sprintf('%.1f ', data.ground_truth.addressed.event_times));
fprintf('  Heard (overheard): %d events at t=[%s]s\n', ...
    n_overheard, sprintf('%.1f ', data.ground_truth.overheard.event_times));
fprintf('  Produced: %d events at t=[%s]s\n', ...
    n_produced, sprintf('%.1f ', data.ground_truth.produced.event_times));
fprintf('  Total events: %d\n', length(data.events));
fprintf('\n');

%% section 2: run pipeline
fprintf('RUNNING GLM PIPELINE...\n');
fprintf('------------------------------------------------------------------------\n');

% bin
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
fprintf('  Binned high gamma: %d bins\n', length(hg_binned.power));

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fprintf('  Design matrix: [%d × %d], %.1f%% sparse\n', ...
    size(X, 1), size(X, 2), 100 * (1 - nnz(X) / numel(X)));

% fit model
fprintf('  Fitting model (λ=%.3f)...\n', data.cfg.model.lambda);
fit = fit_glm_map_hg(X, y, data.cfg);
fprintf('    Converged: %s (%d iterations)\n', mat2str(fit.converged), fit.output.iterations);
fprintf('    Gradient norm: %.6e\n', fit.grad_norm);

% compute metrics
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
fprintf('    R²: %.4f\n', r_squared);

% extract kernels
kernels = extract_kernels(fit.w, colmap, data.stim);
fprintf('  Extracted %d kernel types\n', ...
    sum(structfun(@(x) isstruct(x) && isfield(x, 'weights'), kernels)));
fprintf('\n');

%% section 3: compare kernels
fprintf('KERNEL COMPARISON:\n');
fprintf('------------------------------------------------------------------------\n');

% prepare comparison
kernel_names = {'heard_addressed', 'heard_overheard', 'produced_spontaneous'};

% map to ground truth field names
gt_map = struct();
gt_map.heard_addressed = 'addressed';
gt_map.heard_overheard = 'overheard';
gt_map.produced_spontaneous = 'produced';

% build true_kernels struct for compare_kernels
true_kernels = struct();
for i = 1:length(kernel_names)
    kname = kernel_names{i};
    gt_name = gt_map.(kname);
    true_kernels.(kname) = data.ground_truth.(gt_name);
end

% compare
metrics = compare_kernels(kernels, true_kernels, kernel_names);

%% section 4: visualization
fprintf('CREATING VISUALIZATION...\n');
fprintf('------------------------------------------------------------------------\n');

figure('Position', [100, 100, 1400, 900]);

% subplot 1: full power trace with events
subplot(2, 2, 1);
hold on;

% plot traces
plot(data.stim.t, y, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'DisplayName', 'Observed');
plot(data.stim.t, fit.mu, 'k-', 'LineWidth', 2, 'DisplayName', 'Predicted');

% mark events with different colors
yl = ylim;
for i = 1:length(data.ground_truth.addressed.event_times)
    t_event = data.ground_truth.addressed.event_times(i);
    plot([t_event, t_event], yl, 'b-', 'LineWidth', 1.5);
end
for i = 1:length(data.ground_truth.overheard.event_times)
    t_event = data.ground_truth.overheard.event_times(i);
    plot([t_event, t_event], yl, 'r-', 'LineWidth', 1.5);
end
for i = 1:length(data.ground_truth.produced.event_times)
    t_event = data.ground_truth.produced.event_times(i);
    plot([t_event, t_event], yl, 'g-', 'LineWidth', 1.5);
end

% add dummy plots for legend
plot(NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'Addressed');
plot(NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName', 'Overheard');
plot(NaN, NaN, 'g-', 'LineWidth', 2, 'DisplayName', 'Produced');

title(sprintf('Power Trace with Events (R^2 = %.3f)', r_squared));
xlabel('Time (s)');
ylabel('Power');
legend('Location', 'best');
grid on;
hold off;

% subplot 2: all three kernels overlaid
subplot(2, 2, 2);
hold on;

% define colors
colors = lines(3);

% plot addressed
addr_times = kernels.heard_addressed.times;
addr_fitted = kernels.heard_addressed.weights;
addr_true = data.ground_truth.addressed.kernel_response;

plot(addr_times, addr_fitted, '-', 'Color', colors(1,:), 'LineWidth', 2, ...
    'DisplayName', 'Addressed (fit)');
plot(addr_times, addr_true, '--', 'Color', colors(1,:), 'LineWidth', 1.5, ...
    'DisplayName', 'Addressed (true)');

% plot overheard
over_times = kernels.heard_overheard.times;
over_fitted = kernels.heard_overheard.weights;
over_true = data.ground_truth.overheard.kernel_response;

plot(over_times, over_fitted, '-', 'Color', colors(2,:), 'LineWidth', 2, ...
    'DisplayName', 'Overheard (fit)');
plot(over_times, over_true, '--', 'Color', colors(2,:), 'LineWidth', 1.5, ...
    'DisplayName', 'Overheard (true)');

% plot produced
prod_times = kernels.produced_spontaneous.times;
prod_fitted = kernels.produced_spontaneous.weights;
prod_true = data.ground_truth.produced.kernel_response;

plot(prod_times, prod_fitted, '-', 'Color', colors(3,:), 'LineWidth', 2, ...
    'DisplayName', 'Produced (fit)');
plot(prod_times, prod_true, '--', 'Color', colors(3,:), 'LineWidth', 1.5, ...
    'DisplayName', 'Produced (true)');

title('Kernel Recovery (Fitted vs True)');
xlabel('Time from event (s)');
ylabel('Response amplitude');
legend('Location', 'best');
grid on;
hold off;

% subplot 3: correlation bar plot
subplot(2, 2, 3);

correlations = [...
    metrics.heard_addressed.correlation, ...
    metrics.heard_overheard.correlation, ...
    metrics.produced_spontaneous.correlation];

bar_colors = colors;
b = bar(correlations);
b.FaceColor = 'flat';
b.CData = bar_colors;

hold on;
plot([0, 4], [0.7, 0.7], 'r--', 'LineWidth', 2, 'DisplayName', 'Threshold (0.7)');
hold off;

set(gca, 'XTickLabel', {'Addressed', 'Overheard', 'Produced'});
ylabel('Correlation');
title('Kernel Recovery Quality');
ylim([0, 1]);
grid on;

% add value labels on bars
for i = 1:length(correlations)
    text(i, correlations(i) + 0.05, sprintf('%.3f', correlations(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% subplot 4: residuals
subplot(2, 2, 4);

residuals = y - fit.mu;

yyaxis left;
plot(data.stim.t, residuals, 'k-', 'LineWidth', 0.5);
ylabel('Residual');
hold on;
plot(data.stim.t([1 end]), [0 0], 'r--', 'LineWidth', 1.5);
hold off;

yyaxis right;
histogram(residuals, 30, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
ylabel('Count');

title(sprintf('Residuals (mean=%.3f, std=%.3f)', mean(residuals), std(residuals)));
xlabel('Time (s)');
grid on;

fprintf('Visualization created.\n\n');

%% section 5: checkpoint assessment
fprintf('========================================================================\n');
fprintf('                      CHECKPOINT 2 ASSESSMENT                          \n');
fprintf('========================================================================\n');
fprintf('\n');

% evaluate criteria
all_pass = true;

fprintf('Convergence:\n');
if fit.converged
    fprintf('  [✓] Converged in %d iterations\n', fit.output.iterations);
else
    fprintf('  [✗] Failed to converge\n');
    all_pass = false;
end
fprintf('\n');

fprintf('Overall Fit Quality:\n');
if r_squared > 0.6
    fprintf('  [✓] R²: %.4f (> 0.6)\n', r_squared);
else
    fprintf('  [✗] R²: %.4f (<= 0.6)\n', r_squared);
    all_pass = false;
end
fprintf('\n');

fprintf('Individual Kernel Recovery:\n');
for i = 1:length(kernel_names)
    kname = kernel_names{i};
    if isfield(metrics, kname)
        m = metrics.(kname);

        fprintf('  %s:\n', kname);

        if m.correlation > 0.7
            fprintf('    [✓] Correlation: %.4f (> 0.7)\n', m.correlation);
        else
            fprintf('    [✗] Correlation: %.4f (<= 0.7)\n', m.correlation);
            all_pass = false;
        end

        if m.magnitude_error < 0.4
            fprintf('    [✓] Magnitude error: %.1f%% (< 40%%)\n', m.magnitude_error * 100);
        else
            fprintf('    [✗] Magnitude error: %.1f%% (>= 40%%)\n', m.magnitude_error * 100);
            all_pass = false;
        end
    end
end
fprintf('\n');

fprintf('========================================================================\n');
if all_pass
    fprintf('                                                                        \n');
    fprintf('       ✓ ✓ ✓ CHECKPOINT 2: Multi-event kernel recovery successful! ✓ ✓ ✓\n');
    fprintf('                                                                        \n');
    fprintf('    The model successfully recovers multiple distinct kernels:\n');
    fprintf('      • Heard (addressed): fast positive response\n');
    fprintf('      • Heard (overheard): slow positive response\n');
    fprintf('      • Produced: biphasic response\n');
    fprintf('                                                                        \n');
    fprintf('    All kernels recovered without interference.\n');
    fprintf('    Ready for more complex scenarios (states, basis functions, etc.).\n');
    fprintf('                                                                        \n');
else
    fprintf('                                                                        \n');
    fprintf('    [!] Some criteria not met. Review diagnostics above.\n');
    fprintf('                                                                        \n');
end
fprintf('========================================================================\n');
fprintf('\n');
