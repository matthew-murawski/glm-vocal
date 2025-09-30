function plot_kernels(kernels, outdir)
% section setup
% render kernel weight traces and scalar summaries to disk for quick qc.
if nargin < 2 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fields = fieldnames(kernels);
plotNames = {};
plotData = {};
plotMeta = {};

for ii = 1:numel(fields)
    name = fields{ii};
    value = kernels.(name);
    if isstruct(value) && isfield(value, 'weights')
        plotNames{end+1} = name; %#ok<AGROW>
        plotData{end+1} = value.weights(:); %#ok<AGROW>
        if isfield(value, 'lag_times_s')
            plotMeta{end+1} = value.lag_times_s(:); %#ok<AGROW>
        else
            plotMeta{end+1} = (0:numel(value.weights)-1)'; %#ok<AGROW>
        end
    elseif strcmp(name, 'states') && isstruct(value) && isfield(value, 'weights')
        plotNames{end+1} = name;
        plotData{end+1} = value.weights(:);
        plotMeta{end+1} = [];
    elseif strcmp(name, 'intercept')
        plotNames{end+1} = name;
        plotData{end+1} = kernels.intercept(:);
        plotMeta{end+1} = [];
    end
end

if isempty(plotNames)
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 640, 900]);
try
    nPlots = numel(plotNames);
    if nPlots == 0
        close(fig);
        return
    end

    nCols = min(2, nPlots);
    nRows = ceil(nPlots / nCols);
    tile = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for jj = 1:nPlots
        ax = nexttile(tile);
        data = plotData{jj};
        meta = plotMeta{jj};
        name = plotNames{jj};
        if strcmp(name, 'states')
            bar(ax, data);
            xticks(ax, 1:numel(data));
            if isfield(kernels.states, 'names')
                xticklabels(ax, kernels.states.names);
            else
                xticklabels(ax, compose('state %d', 1:numel(data)));
            end
            ylabel(ax, 'weight');
        elseif strcmp(name, 'intercept')
            yline(ax, data, 'LineWidth', 1.5);
            xlim(ax, [-0.5, 0.5]);
            ylabel(ax, 'intercept');
            xticklabels(ax, {});
            pad = max(0.05, 0.1 * max(1, abs(data)));
            ylim(ax, data + [-pad, pad]);
        else
            plot(ax, meta, data, 'LineWidth', 1.5);
            xlabel(ax, 'lag (s)');
            ylabel(ax, 'weight');
        end
        title(ax, strrep(name, '_', ' '));
        grid(ax, 'on');
        if ~strcmp(name, 'intercept')
            axis(ax, 'tight');
        end
        set(ax, 'PlotBoxAspectRatio', [1 1 1]);
    end
    filepath = fullfile(outdir, 'kernels.pdf');
    exportgraphics(fig, filepath, 'ContentType', 'vector');
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end
