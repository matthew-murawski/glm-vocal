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

fig = figure('Visible', 'off');
try
    nPlots = numel(plotNames);
    tiledlayout(fig, nPlots, 1, 'TileSpacing', 'compact');
    for jj = 1:nPlots
        nexttile;
        data = plotData{jj};
        meta = plotMeta{jj};
        name = plotNames{jj};
        if strcmp(name, 'states')
            bar(data);
            xticks(1:numel(data));
            if isfield(kernels.states, 'names')
                xticklabels(kernels.states.names);
            else
                xticklabels(compose('state %d', 1:numel(data)));
            end
            ylabel('weight');
        elseif strcmp(name, 'intercept')
            plot(0, data, 'o');
            xlim([-0.5, 0.5]);
            ylabel('intercept');
            xticklabels({});
        else
            plot(meta, data, '-o');
            xlabel('lag (s)');
            ylabel('weight');
        end
        title(strrep(name, '_', ' '));
        grid on;
    end
    filepath = fullfile(outdir, 'kernels.png');
    saveas(fig, filepath);
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end
