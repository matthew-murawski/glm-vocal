function fig = plot_lfp_traces(t, lfp_actual, lfp_pred, channel_id, metricsOut)
%PLOT_LFP_TRACES Plot actual vs predicted LFP traces.
%   fig = PLOT_LFP_TRACES(t, lfp_actual, lfp_pred, channel_id, metricsOut)
%   creates a figure showing the actual and predicted LFP signals overlaid,
%   along with evaluation metrics.
%
%   Inputs:
%     t           - time vector (n_bins × 1) in seconds
%     lfp_actual  - actual LFP signal (n_bins × 1)
%     lfp_pred    - predicted LFP signal (n_bins × 1)
%     channel_id  - (optional) channel identifier for title
%     metricsOut  - (optional) struct from metrics_lfp with .r2, .correlation, etc.
%
%   Returns:
%     fig - figure handle

% section input validation
if nargin < 4 || isempty(channel_id)
    channel_id = 'Unknown';
end
if nargin < 5
    metricsOut = struct();
end

t = t(:);
lfp_actual = lfp_actual(:);
lfp_pred = lfp_pred(:);

% section figure creation
fig = figure('Position', [100, 100, 1200, 400]);

% plot actual vs predicted traces
plot(t, lfp_actual, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Actual LFP');
hold on;
plot(t, lfp_pred, 'r-', 'LineWidth', 1.0, 'DisplayName', 'Predicted LFP');
hold off;

xlabel('Time (s)', 'FontSize', 12);
ylabel('LFP (normalized)', 'FontSize', 12);

% construct title with metrics
title_str = sprintf('LFP Prediction: Channel %s', num2str(channel_id));
if isfield(metricsOut, 'r2') && ~isnan(metricsOut.r2)
    title_str = sprintf('%s | R² = %.3f', title_str, metricsOut.r2);
end
if isfield(metricsOut, 'correlation') && ~isnan(metricsOut.correlation)
    title_str = sprintf('%s | Corr = %.3f', title_str, metricsOut.correlation);
end
if isfield(metricsOut, 'rmse') && ~isnan(metricsOut.rmse)
    title_str = sprintf('%s | RMSE = %.3f', title_str, metricsOut.rmse);
end
title(title_str, 'FontSize', 14, 'Interpreter', 'none');

legend('Location', 'best', 'FontSize', 10);
grid on;

% set axis limits
xlim([min(t), max(t)]);

end
