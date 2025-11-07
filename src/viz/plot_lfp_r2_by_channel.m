function fig = plot_lfp_r2_by_channel(results)
%PLOT_LFP_R2_BY_CHANNEL Create a bar plot of R² values across channels.
%   fig = PLOT_LFP_R2_BY_CHANNEL(results) creates a bar plot showing
%   the R-squared value for each LFP channel, allowing comparison of
%   model fit quality across electrodes.
%
%   Inputs:
%     results - array of structs (one per channel) with fields:
%               .channel_id - channel identifier
%               .metrics    - struct from metrics_lfp with .r2 field
%
%   Returns:
%     fig - figure handle

% section input validation
if ~isstruct(results) || isempty(results)
    error('plot_lfp_r2_by_channel:InvalidInput', 'results must be a non-empty struct array.');
end

n_channels = numel(results);

% section extract R² values and channel IDs
r2_values = zeros(n_channels, 1);
channel_ids = cell(n_channels, 1);

for ch = 1:n_channels
    if isfield(results(ch), 'metrics') && isfield(results(ch).metrics, 'r2')
        r2_values(ch) = results(ch).metrics.r2;
    else
        r2_values(ch) = NaN;
    end

    if isfield(results(ch), 'channel_id')
        channel_ids{ch} = num2str(results(ch).channel_id);
    else
        channel_ids{ch} = num2str(ch);
    end
end

% section figure creation
fig = figure('Position', [100, 100, 1200, 600]);

% create bar plot
bar(1:n_channels, r2_values, 'FaceColor', [0.2, 0.4, 0.8]);
hold on;

% add horizontal line at R² = 0 for reference
plot([0.5, n_channels + 0.5], [0, 0], 'k--', 'LineWidth', 1);

hold off;

% labels and title
xlabel('Channel', 'FontSize', 12);
ylabel('R² (Fraction of Variance Explained)', 'FontSize', 12);
title('LFP Prediction Performance Across Channels', 'FontSize', 14);

% set x-axis tick labels to channel IDs
set(gca, 'XTick', 1:n_channels, 'XTickLabel', channel_ids);
if n_channels > 20
    % rotate labels if many channels
    set(gca, 'XTickLabelRotation', 45);
end

% set y-axis limits
ylim_min = min(-0.1, min(r2_values) - 0.05);
ylim_max = max(1.0, max(r2_values) + 0.05);
ylim([ylim_min, ylim_max]);

grid on;

% compute summary statistics
mean_r2 = nanmean(r2_values);
median_r2 = nanmedian(r2_values);

% add text annotation with summary stats
annotation_str = sprintf('Mean R² = %.3f | Median R² = %.3f', mean_r2, median_r2);
text(0.5, 0.98, annotation_str, 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'BackgroundColor', 'white', 'EdgeColor', 'black');

end
