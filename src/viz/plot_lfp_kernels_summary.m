function fig = plot_lfp_kernels_summary(kernels, channel_id, metrics)
%PLOT_LFP_KERNELS_SUMMARY Create summary plot of all fitted kernels for a channel.
%   fig = PLOT_LFP_KERNELS_SUMMARY(kernels, channel_id, metrics) creates
%   a multi-panel figure showing the shape of all fitted predictive kernels
%   for a single LFP channel.
%
%   Inputs:
%     kernels    - struct with kernel fields (from unpack_params)
%     channel_id - identifier for this channel
%     metrics    - struct with fit metrics (e.g., r2, correlation)
%
%   Returns:
%     fig - figure handle

% create figure
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% determine which kernels to plot
kernel_names = {};
kernel_data = {};
kernel_times = {};
kernel_colors = {};

% color scheme
color_produced = [0.8, 0.2, 0.2];      % red for produced
color_addressed = [0.2, 0.6, 0.8];     % blue for addressed
color_overheard = [0.4, 0.7, 0.4];     % green for overheard
color_history = [0.6, 0.4, 0.8];       % purple for history
color_default = [0.3, 0.3, 0.3];       % gray for others

% collect produced kernels
if isfield(kernels, 'produced_fields')
    for ii = 1:numel(kernels.produced_fields)
        fname = kernels.produced_fields{ii};
        if isfield(kernels, fname) && isfield(kernels.(fname), 'weights')
            kernel_names{end+1} = strrep(fname, '_', ' '); %#ok<AGROW>
            kernel_data{end+1} = kernels.(fname).weights(:); %#ok<AGROW>
            if isfield(kernels.(fname), 'lag_times_s')
                kernel_times{end+1} = kernels.(fname).lag_times_s(:); %#ok<AGROW>
            else
                kernel_times{end+1} = (0:numel(kernels.(fname).weights)-1)'; %#ok<AGROW>
            end
            kernel_colors{end+1} = color_produced; %#ok<AGROW>
        end
    end
end

% collect heard kernels
if isfield(kernels, 'heard_fields')
    for ii = 1:numel(kernels.heard_fields)
        fname = kernels.heard_fields{ii};
        if isfield(kernels, fname) && isfield(kernels.(fname), 'weights')
            kernel_names{end+1} = strrep(fname, '_', ' '); %#ok<AGROW>
            kernel_data{end+1} = kernels.(fname).weights(:); %#ok<AGROW>
            if isfield(kernels.(fname), 'lag_times_s')
                kernel_times{end+1} = kernels.(fname).lag_times_s(:); %#ok<AGROW>
            else
                kernel_times{end+1} = (0:numel(kernels.(fname).weights)-1)'; %#ok<AGROW>
            end

            % assign color based on type
            if contains(fname, 'addressed')
                kernel_colors{end+1} = color_addressed; %#ok<AGROW>
            elseif contains(fname, 'overheard')
                kernel_colors{end+1} = color_overheard; %#ok<AGROW>
            else
                kernel_colors{end+1} = color_default; %#ok<AGROW>
            end
        end
    end
end

% collect spike/LFP history kernel if present
if isfield(kernels, 'spike_history') && isfield(kernels.spike_history, 'weights')
    kernel_names{end+1} = 'history'; %#ok<AGROW>
    kernel_data{end+1} = kernels.spike_history.weights(:); %#ok<AGROW>
    if isfield(kernels.spike_history, 'lag_times_s')
        kernel_times{end+1} = kernels.spike_history.lag_times_s(:); %#ok<AGROW>
    else
        kernel_times{end+1} = (0:numel(kernels.spike_history.weights)-1)'; %#ok<AGROW>
    end
    kernel_colors{end+1} = color_history; %#ok<AGROW>
end

% determine layout
n_kernels = numel(kernel_names);
if n_kernels == 0
    close(fig);
    warning('No kernels found to plot');
    return;
end

n_cols = min(3, n_kernels);
n_rows = ceil(n_kernels / n_cols);

% create tiled layout
tile = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

% plot each kernel
for ii = 1:n_kernels
    ax = nexttile(tile);

    times = kernel_times{ii};
    weights = kernel_data{ii};
    color = kernel_colors{ii};

    % plot zero line
    hold(ax, 'on');
    plot(ax, times, zeros(size(times)), ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);

    % plot kernel
    plot(ax, times, weights, 'LineWidth', 2.5, 'Color', color);

    % add vertical line at t=0 for event-locked kernels
    if min(times) < 0 && max(times) > 0
        xline(ax, 0, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5);
    end

    hold(ax, 'off');

    % formatting
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Weight');
    title(ax, kernel_names{ii}, 'Interpreter', 'none');
    grid(ax, 'on');
    axis(ax, 'tight');

    % add subtle shading for positive/negative regions
    ylims = ylim(ax);
    xlims = xlim(ax);
    ylim(ax, ylims);
    xlim(ax, xlims);
end

% add overall title with channel info and metrics
title_str = sprintf('Channel %s - Fitted Kernels', num2str(channel_id));
if ~isempty(metrics)
    if isfield(metrics, 'r2')
        title_str = sprintf('%s (R² = %.3f)', title_str, metrics.r2);
    end
    if isfield(metrics, 'correlation')
        title_str = sprintf('%s, Corr = %.3f', title_str, metrics.correlation);
    end
end
title(tile, title_str, 'FontSize', 14, 'FontWeight', 'bold');

end
