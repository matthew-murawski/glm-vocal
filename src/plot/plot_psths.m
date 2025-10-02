function plot_psths(sp, heardEvents, producedEvents, cfg, outdir)
% render psths for heard calls plus produced calls split by conversational context.

% section setup figure
% build a four-panel layout so heard responses and produced subtypes are visible together.
combinedEvents = [heardEvents(:); producedEvents(:)];
if ~isempty(combinedEvents)
    [~, order] = sort(double([combinedEvents.t_on]));
    combinedEvents = combinedEvents(order);
end

if isempty(combinedEvents)
    producedMask = false(0, 1);
else
    kindCells = cellfun(@char, {combinedEvents.kind}, 'UniformOutput', false);
    producedMask = strcmp(kindCells, 'produced');
end
lookbackWindowS = 5.0;
categoryIdx = classify_produced_events(combinedEvents, producedMask, lookbackWindowS, cfg);
producedNames = categoryIdx.names;
nProduced = numel(producedNames);

totalPanels = 1 + max(nProduced, 1);
nCols = min(2, totalPanels);
nRows = ceil(totalPanels / nCols);

fig = figure('Name', 'PSTHs', 'Color', 'w', 'Position', [100, 100, 900, 700], 'Visible', 'off');
tile = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact');

% section heard psth
% first panel shows alignment to heard call onsets.
ax = nexttile(tile);
if ~isempty(heardEvents)
    psth(sp.spike_times, [heardEvents.t_on], cfg.heard_window_s, cfg.dt, ...
        'Axes', ax, 'Title', 'psth: heard calls');
else
    title(ax, 'psth: heard calls (no events)');
    box(ax, 'on');
    set(ax, 'xtick', [], 'ytick', []);
end

% section produced psths
% remaining panels show each produced-call context separately.
for kk = 1:nProduced
    ax = nexttile(tile);
    fieldName = producedNames{kk};
    idx = categoryIdx.(fieldName);
    panelTitle = sprintf('psth: %s', strrep(fieldName, '_', ' '));
    if isempty(idx)
        title(ax, [panelTitle, ' (no events)']);
        box(ax, 'on');
        set(ax, 'xtick', [], 'ytick', []);
        continue
    end

    producedSubset = combinedEvents(idx);
    psth(sp.spike_times, [producedSubset.t_on], cfg.produced_window_s, cfg.dt, ...
        'Axes', ax, 'Title', panelTitle);
end

% section save figure
% we save the whole figure as a pdf in the specified output directory.
output_path = fullfile(outdir, 'psth.pdf');
try
    % using print to save as a high-quality pdf.
    print(fig, output_path, '-dpdf', '-r300', '-bestfit');
catch ME
    warning('could not save psth plot to %s. error: %s', output_path, ME.message);
    % close the figure if saving failed to avoid it popping up.
    close(fig);
    rethrow(ME);
end

% close the figure now that we're done with it.
close(fig);

end
