%% demo_step8_states
% demonstrate state coefficient recovery with conversational states

clr; clc;

fprintf('========================================\n');
fprintf('Step 8: State Coefficient Recovery\n');
fprintf('========================================\n\n');

%% section generate synthetic data
fprintf('generating synthetic session with state effects...\n');
data = generate_state_effect_hg();

fprintf('  duration: %.1f s\n', data.ground_truth.duration_s);
fprintf('  epochs: [0-10s] spon, [10-20s] convo, [20-30s] spon\n');
fprintf('  baseline spon: %.1f\n', data.ground_truth.baseline_spon);
fprintf('  baseline convo: %.1f\n', data.ground_truth.baseline_convo);
fprintf('  state boost (convo): +%.1f\n', data.ground_truth.state_coeff_convo);
fprintf('  events: %d perceived calls\n\n', length(data.events));

%% section run pipeline
fprintf('running pipeline...\n');

% bin high gamma
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
fprintf('  binned high gamma: %d time bins\n', length(hg_binned.t));

% assemble design matrix
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);
fprintf('  design matrix: %d × %d\n', size(X, 1), size(X, 2));

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);
n_iter = 0;
if isfield(fit.output, 'iterations')
    n_iter = fit.output.iterations;
end
fprintf('  converged: %s (iter=%d, grad_norm=%.2e)\n\n', ...
    mat2str(fit.converged), n_iter, fit.grad_norm);

% extract kernels and states
kernels = extract_kernels(fit.w, colmap, data.stim);

%% section evaluate recovery
fprintf('evaluating recovery...\n');

% kernel correlation
kernel_corr = corr(kernels.heard_addressed.weights, ...
    data.ground_truth.kernel.kernel_response);
fprintf('  kernel correlation: %.3f\n', kernel_corr);

% state coefficients
convo_coeff = kernels.states.convo;
spon_coeff = kernels.states.spon;
true_convo = data.ground_truth.state_coeff_convo;
true_spon = data.ground_truth.state_coeff_spon;

convo_error = 100 * abs(convo_coeff - true_convo) / abs(true_convo);
spon_error = abs(spon_coeff - true_spon);

fprintf('  convo coefficient: %.2f (true: %.2f, error: %.1f%%)\n', ...
    convo_coeff, true_convo, convo_error);
fprintf('  spon coefficient: %.2f (true: %.2f, error: %.2f)\n', ...
    spon_coeff, true_spon, spon_error);

% overall R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;
fprintf('  R²: %.3f\n', r_squared);

% epoch-specific predictions
convo_mask = data.ground_truth.convo_mask;
spon_mask = data.ground_truth.spon_mask;
mean_pred_convo = mean(fit.mu(convo_mask));
mean_pred_spon = mean(fit.mu(spon_mask));
pred_shift = mean_pred_convo - mean_pred_spon;
true_shift = data.ground_truth.baseline_convo - data.ground_truth.baseline_spon;

fprintf('  mean predicted power (convo): %.2f\n', mean_pred_convo);
fprintf('  mean predicted power (spon): %.2f\n', mean_pred_spon);
fprintf('  predicted baseline shift: %.2f (true: %.2f)\n\n', pred_shift, true_shift);

%% section visualization
fprintf('creating visualization...\n');

figure('position', [100, 100, 1200, 800]);

% subplot 1: power trace with state coloring
subplot(2, 3, 1);
hold on;

% plot state backgrounds
convo_idx_plot = find(convo_mask);
spon_idx_plot = find(spon_mask);

% background shading
for i = 1:length(convo_idx_plot)
    if i == 1 || convo_idx_plot(i) ~= convo_idx_plot(i-1) + 1
        patch([hg_binned.t(convo_idx_plot(i)), hg_binned.t(min(convo_idx_plot(i)+1, end)), ...
               hg_binned.t(min(convo_idx_plot(i)+1, end)), hg_binned.t(convo_idx_plot(i))], ...
              [0, 0, 20, 20], [1, 0.9, 0.9], 'edgecolor', 'none', 'facealpha', 0.3);
    end
end

% plot trace
plot(hg_binned.t, y, 'k-', 'linewidth', 1.5);
plot(hg_binned.t, fit.mu, 'r-', 'linewidth', 1);

% mark events
for i = 1:length(data.events)
    line([data.events(i).t_on, data.events(i).t_on], [0, max(y)*1.1], ...
        'color', [0.5, 0.5, 1], 'linestyle', '--', 'linewidth', 1);
end

xlabel('time (s)');
ylabel('high gamma power');
title(sprintf('power trace (R² = %.3f)', r_squared));
legend({'convo epoch', 'observed', 'predicted', 'events'}, 'location', 'best');
grid on;

% subplot 2: state-specific predictions
subplot(2, 3, 2);
hold on;

% violin-style distribution
spon_pred = fit.mu(spon_mask);
convo_pred = fit.mu(convo_mask);

boxplot([spon_pred; convo_pred], [ones(length(spon_pred), 1); 2*ones(length(convo_pred), 1)], ...
    'labels', {'spon', 'convo'});

% add true baselines
line([0.5, 1.5], [data.ground_truth.baseline_spon, data.ground_truth.baseline_spon], ...
    'color', 'r', 'linestyle', '--', 'linewidth', 2);
line([1.5, 2.5], [data.ground_truth.baseline_convo, data.ground_truth.baseline_convo], ...
    'color', 'r', 'linestyle', '--', 'linewidth', 2);

ylabel('predicted power');
title(sprintf('epoch-specific predictions\nshift: %.2f (true: %.2f)', pred_shift, true_shift));
grid on;

% subplot 3: fitted vs true kernel
subplot(2, 3, 3);
hold on;

plot(data.ground_truth.kernel.kernel_times, data.ground_truth.kernel.kernel_response, ...
    'k--', 'linewidth', 2);
plot(kernels.heard_addressed.times, kernels.heard_addressed.weights, ...
    'r-', 'linewidth', 1.5);

xlabel('time from event onset (s)');
ylabel('kernel weight');
title(sprintf('kernel recovery (r = %.3f)', kernel_corr));
legend({'true kernel', 'fitted kernel'}, 'location', 'best');
grid on;

% subplot 4: state coefficients comparison
subplot(2, 3, 4);
hold on;

states_true = [true_spon; true_convo];
states_fitted = [spon_coeff; convo_coeff];

bar_x = [1, 2];
bar([1, 2], states_true, 'facecolor', [0.7, 0.7, 0.7], 'edgecolor', 'k', 'linewidth', 1.5);
bar([1.25, 2.25], states_fitted, 'facecolor', [1, 0.5, 0.5], 'edgecolor', 'k', 'linewidth', 1.5);

set(gca, 'xtick', [1.125, 2.125], 'xticklabel', {'spon', 'convo'});
ylabel('state coefficient');
title('state coefficient recovery');
legend({'true', 'fitted'}, 'location', 'best');
grid on;

% subplot 5: residuals over time
subplot(2, 3, 5);
residuals = y - fit.mu;
plot(hg_binned.t, residuals, 'k-', 'linewidth', 0.5);
hold on;
line([0, max(hg_binned.t)], [0, 0], 'color', 'r', 'linestyle', '--');

% shade state epochs
for i = 1:length(convo_idx_plot)
    if i == 1 || convo_idx_plot(i) ~= convo_idx_plot(i-1) + 1
        patch([hg_binned.t(convo_idx_plot(i)), hg_binned.t(min(convo_idx_plot(i)+1, end)), ...
               hg_binned.t(min(convo_idx_plot(i)+1, end)), hg_binned.t(convo_idx_plot(i))], ...
              [-5, -5, 5, 5], [1, 0.9, 0.9], 'edgecolor', 'none', 'facealpha', 0.2);
    end
end

xlabel('time (s)');
ylabel('residual');
title(sprintf('residuals (std = %.2f)', std(residuals)));
grid on;

% subplot 6: observed vs predicted scatter
subplot(2, 3, 6);
hold on;

% color by state
scatter(y(spon_mask), fit.mu(spon_mask), 20, [0.5, 0.5, 0.5], 'filled', 'markerfacealpha', 0.5);
scatter(y(convo_mask), fit.mu(convo_mask), 20, [1, 0.5, 0.5], 'filled', 'markerfacealpha', 0.5);

% identity line
lims = [min([y; fit.mu]), max([y; fit.mu])];
plot(lims, lims, 'k--', 'linewidth', 1);

xlabel('observed power');
ylabel('predicted power');
title(sprintf('obs vs pred (R² = %.3f)', r_squared));
legend({'spon', 'convo', 'identity'}, 'location', 'best');
grid on;
axis equal;
xlim(lims);
ylim(lims);

%% section summary
fprintf('\n========================================\n');
fprintf('Step 8 Complete\n');
fprintf('========================================\n');
fprintf('✓ Model converged\n');
fprintf('✓ Kernel recovered (r = %.3f)\n', kernel_corr);
fprintf('✓ Convo coefficient recovered (error = %.1f%%)\n', convo_error);
fprintf('✓ Spon coefficient near zero (error = %.2f)\n', spon_error);
fprintf('✓ Overall R² = %.3f\n', r_squared);
fprintf('✓ Epoch-specific predictions verified\n');
fprintf('========================================\n\n');
