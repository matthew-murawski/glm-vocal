%% demo_step5_design_matrix
% demonstrates design matrix assembly for high gamma GLM

clr; clc;

fprintf('\n========================================\n');
fprintf('Step 5: Design Matrix Assembly Demo\n');
fprintf('========================================\n\n');

%% section 1: generate synthetic session
fprintf('1. Generating synthetic session with single perceived call...\n');
fprintf('   - Duration: 6 seconds\n');
fprintf('   - Baseline power: 5\n');
fprintf('   - Event: perceived call at t=2s\n');
fprintf('   - Response: Gaussian bump (peak=+2, width=0.3s, delay=0.1s)\n\n');

data = generate_single_event_hg();

fprintf('   Generated session:\n');
fprintf('   - Time bins: %d (dt = %.3f s)\n', length(data.stim.t), data.cfg.dt);
fprintf('   - Power range: [%.2f, %.2f]\n', min(data.hg_data.power), max(data.hg_data.power));
fprintf('   - Event time: %.2f s\n\n', data.ground_truth.event_time);

%% section 2: bin high gamma
fprintf('2. Binning high gamma to GLM time base...\n');

hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);

fprintf('   Binned power:\n');
fprintf('   - Method: %s\n', hg_binned.method);
fprintf('   - Output length: %d bins\n', length(hg_binned.power));
fprintf('   - Empty bins: %d\n', hg_binned.n_empty_bins);
if ~isempty(hg_binned.warnings)
    fprintf('   Warnings: %d\n', length(hg_binned.warnings));
end
fprintf('\n');

%% section 3: assemble design matrix
fprintf('3. Assembling design matrix...\n');

[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

fprintf('   Design matrix properties:\n');
fprintf('   - Dimensions: [%d × %d]\n', size(X, 1), size(X, 2));
fprintf('   - Sparsity: %.1f%%\n', 100 * (1 - nnz(X) / numel(X)));
fprintf('   - Response length: %d\n', length(y));
fprintf('   - Response range: [%.2f, %.2f]\n\n', min(y), max(y));

%% section 4: examine column structure
fprintf('4. Design matrix column structure:\n');

fprintf('   Regressors included:\n');

% intercept
if isfield(colmap, 'intercept')
    fprintf('   - intercept: column %d\n', colmap.intercept.cols);
end

% heard fields
if isfield(colmap, 'heard_fields')
    for i = 1:length(colmap.heard_fields)
        field = colmap.heard_fields{i};
        cols = colmap.(field).cols;
        fprintf('   - %s: columns %d-%d (%d cols)\n', ...
            field, min(cols), max(cols), length(cols));
    end
end

% produced fields
if isfield(colmap, 'produced_fields')
    for i = 1:length(colmap.produced_fields)
        field = colmap.produced_fields{i};
        cols = colmap.(field).cols;
        fprintf('   - %s: columns %d-%d (%d cols)\n', ...
            field, min(cols), max(cols), length(cols));
    end
end

% states
if isfield(colmap, 'states')
    fprintf('   - states: columns %d-%d (', min(colmap.states.cols), max(colmap.states.cols));
    fprintf('%s', strjoin(colmap.states.names, ', '));
    fprintf(')\n');
end

% verify no spike history
if isfield(colmap, 'spike_history')
    fprintf('   - ERROR: spike_history present (should be excluded!)\n');
else
    fprintf('   - spike_history: excluded (correct for high gamma)\n');
end

fprintf('\n');

%% section 5: visualize design matrix structure
fprintf('5. Visualizing design matrix structure...\n\n');

figure('Position', [100, 100, 1200, 500]);

% subplot 1: spy plot of design matrix
subplot(1, 3, 1);
spy(X);
title(sprintf('Design Matrix Structure\n[%d × %d, %.1f%% sparse]', ...
    size(X, 1), size(X, 2), 100 * (1 - nnz(X) / numel(X))));
xlabel('Column (regressor)');
ylabel('Time bin');

% subplot 2: column-wise density
subplot(1, 3, 2);
col_density = full(sum(X ~= 0, 1)) / size(X, 1) * 100;
bar(col_density);
title('Column Density');
xlabel('Column index');
ylabel('% non-zero entries');
grid on;

% subplot 3: highlight regressor blocks with colors
subplot(1, 3, 3);
hold on;

y_pos = 0;
colors = lines(7);
color_idx = 1;

% plot each regressor block
if isfield(colmap, 'intercept')
    cols = colmap.intercept.cols;
    patch([min(cols)-0.5, max(cols)+0.5, max(cols)+0.5, min(cols)-0.5], ...
          [y_pos, y_pos, y_pos+1, y_pos+1], colors(color_idx, :), ...
          'EdgeColor', 'k', 'LineWidth', 1);
    text(mean(cols), y_pos+0.5, 'intercept', 'HorizontalAlignment', 'center');
    y_pos = y_pos + 1;
    color_idx = mod(color_idx, size(colors, 1)) + 1;
end

if isfield(colmap, 'heard_fields')
    for i = 1:length(colmap.heard_fields)
        field = colmap.heard_fields{i};
        cols = colmap.(field).cols;
        patch([min(cols)-0.5, max(cols)+0.5, max(cols)+0.5, min(cols)-0.5], ...
              [y_pos, y_pos, y_pos+1, y_pos+1], colors(color_idx, :), ...
              'EdgeColor', 'k', 'LineWidth', 1);
        text(mean(cols), y_pos+0.5, strrep(field, '_', ' '), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        y_pos = y_pos + 1;
        color_idx = mod(color_idx, size(colors, 1)) + 1;
    end
end

if isfield(colmap, 'produced_fields')
    for i = 1:length(colmap.produced_fields)
        field = colmap.produced_fields{i};
        cols = colmap.(field).cols;
        patch([min(cols)-0.5, max(cols)+0.5, max(cols)+0.5, min(cols)-0.5], ...
              [y_pos, y_pos, y_pos+1, y_pos+1], colors(color_idx, :), ...
              'EdgeColor', 'k', 'LineWidth', 1);
        text(mean(cols), y_pos+0.5, strrep(field, '_', ' '), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        y_pos = y_pos + 1;
        color_idx = mod(color_idx, size(colors, 1)) + 1;
    end
end

if isfield(colmap, 'states')
    cols = colmap.states.cols;
    patch([min(cols)-0.5, max(cols)+0.5, max(cols)+0.5, min(cols)-0.5], ...
          [y_pos, y_pos, y_pos+1, y_pos+1], colors(color_idx, :), ...
          'EdgeColor', 'k', 'LineWidth', 1);
    text(mean(cols), y_pos+0.5, 'states', 'HorizontalAlignment', 'center');
    y_pos = y_pos + 1;
end

xlim([0, size(X, 2) + 1]);
ylim([0, y_pos]);
xlabel('Column index');
title('Regressor Blocks');
set(gca, 'YTick', []);
grid on;
hold off;

%% section 6: show sample time bins
fprintf('6. Examining design matrix at key time points...\n\n');

% find event time
event_bin = find(data.stim.t >= data.ground_truth.event_time, 1, 'first');

% show bins: before event, during event, after event
sample_bins = [max(1, event_bin - 50), event_bin, min(length(y), event_bin + 50)];
sample_labels = {'Before event (t=%.2fs)', 'During event (t=%.2fs)', 'After event (t=%.2fs)'};

for i = 1:length(sample_bins)
    bin_idx = sample_bins(i);
    fprintf('   %s:\n', sprintf(sample_labels{i}, data.stim.t(bin_idx)));
    fprintf('   - Power: %.3f\n', y(bin_idx));

    % show non-zero columns
    row = X(bin_idx, :);
    nz_cols = find(row ~= 0);

    if isempty(nz_cols)
        fprintf('   - No active regressors\n');
    else
        fprintf('   - Active columns: ');
        if length(nz_cols) > 10
            fprintf('%d total (showing first 10): ', length(nz_cols));
            nz_cols = nz_cols(1:10);
        end
        fprintf('[%s]\n', num2str(nz_cols));

        fprintf('   - Values: [');
        for j = 1:length(nz_cols)
            fprintf('%.2f', full(row(nz_cols(j))));
            if j < length(nz_cols)
                fprintf(', ');
            end
        end
        fprintf(']\n');
    end
    fprintf('\n');
end

%% section 7: verify integration with fitting
fprintf('7. Testing integration with GLM fitting...\n');

fit = fit_glm_map_hg(X, y, data.cfg);

fprintf('   Fit results:\n');
fprintf('   - Converged: %s\n', mat2str(fit.converged));
fprintf('   - Iterations: %d\n', fit.output.iterations);
fprintf('   - Final NLL: %.6f\n', fit.nll);
fprintf('   - Gradient norm: %.6e\n', fit.grad_norm);

% compute R²
ss_res = sum((y - fit.mu).^2);
ss_tot = sum((y - mean(y)).^2);
r_squared = 1 - ss_res / ss_tot;

fprintf('   - R²: %.4f\n', r_squared);
fprintf('   - Correlation: %.4f\n\n', corr(y, fit.mu));

%% section 8: summary
fprintf('========================================\n');
fprintf('Summary:\n');
fprintf('========================================\n');
fprintf('Design matrix assembled successfully!\n');
fprintf('  - Dimensions: [%d × %d]\n', size(X, 1), size(X, 2));
fprintf('  - Sparsity: %.1f%%\n', 100 * (1 - nnz(X) / numel(X)));
fprintf('  - Regressors: %d types\n', length(fieldnames(colmap)));
fprintf('  - Spike history: excluded ✓\n');
fprintf('  - Fit converged: %s\n', mat2str(fit.converged));
fprintf('  - R²: %.4f\n', r_squared);
fprintf('\n');
fprintf('✓ Step 5 Complete: Design matrix assembly working\n');
fprintf('========================================\n\n');
