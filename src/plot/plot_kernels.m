function plot_kernels(kernels, ptest, outdir, outfile)
% section setup
% render kernel weight traces and scalar summaries to disk for quick qc.
if nargin < 2
    ptest = [];
    outdir = pwd;
elseif nargin < 3
    if ischar(ptest) || isstring(ptest)
        outdir = char(ptest);
        ptest = [];
    else
        outdir = pwd;
    end
end
if nargin == 3 && (ischar(outdir) || isstring(outdir))
    outdir = char(outdir);
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% choose output filename; default keeps backward compatibility.
if nargin < 4 || isempty(outfile)
    outfile = 'kernels.pdf';
else
    if isstring(outfile) || ischar(outfile)
        outfile = char(outfile);
    else
        error('plot_kernels:InvalidOutfile', 'outfile must be a char or string.');
    end
    % ensure a pdf suffix for consistency
    [~, ~, ext] = fileparts(outfile);
    if isempty(ext)
        outfile = [outfile, '.pdf']; %#ok<AGROW>
    end
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
            lineStyle = '-';
            lineColor = [0, 0, 0];
            ciLower = [];
            ciUpper = [];
            if ~isempty(ptest) && isfield(ptest, name)
                if ptest.(name).p_value >= 0.05
                    lineStyle = '--';
                    lineColor = [0.5, 0.5, 0.5];
                end

                ciLower = ptest.(name).ci_lower(:);
                ciUpper = ptest.(name).ci_upper(:);
            end

            hold(ax, 'on');
            if ~isempty(ciLower) && ~isempty(ciUpper)
                fill(ax, [meta; flipud(meta)], [ciLower; flipud(ciUpper)], ...
                    [0.7 0.7 0.7], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            end
            if ~isempty(meta)
                plot(ax, meta, zeros(size(meta)), ':', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);
            end
            plot(ax, meta, data, 'LineWidth', 1.5, 'LineStyle', lineStyle, 'Color', lineColor);
            hold(ax, 'off');
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
    filepath = fullfile(outdir, outfile);
    exportgraphics(fig, filepath, 'ContentType', 'vector');
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end
