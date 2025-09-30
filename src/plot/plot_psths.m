function plot_psths(sp, heardEvents, producedEvents, cfg, outdir)
% render psths for heard calls plus produced calls split by conversational context.

% section setup figure
% build a four-panel layout so heard responses and produced subtypes are visible together.
fig = figure('Name', 'PSTHs', 'Color', 'w', 'Position', [100, 100, 900, 700], 'Visible', 'off');
tile = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');

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

% section prepare produced subsets
% categorise produced calls using the same five second lookback as the glm preprocessing.
combinedEvents = [heardEvents(:); producedEvents(:)];
if ~isempty(combinedEvents)
    [~, order] = sort(double([combinedEvents.t_on]));
    combinedEvents = combinedEvents(order);
end

if isempty(combinedEvents)
    kindCells = {};
    producedMask = false(0, 1);
else
    kindCells = cellfun(@char, {combinedEvents.kind}, 'UniformOutput', false);
    producedMask = strcmp(kindCells, 'produced');
end
lookbackWindowS = 5.0;
categoryIdx = classify_produced_events(combinedEvents, producedMask, lookbackWindowS);

producedSets = {
    struct('title', 'psth: produced spont', 'indices', categoryIdx.produced_spontaneous), ...
    struct('title', 'psth: produced after heard', 'indices', categoryIdx.produced_after_heard), ...
    struct('title', 'psth: produced after produced', 'indices', categoryIdx.produced_after_produced)
};

% section produced psths
% remaining panels show each produced-call context separately.
for kk = 1:numel(producedSets)
    ax = nexttile(tile);
    idx = producedSets{kk}.indices;
    if isempty(idx)
        title(ax, [producedSets{kk}.title, ' (no events)']);
        box(ax, 'on');
        set(ax, 'xtick', [], 'ytick', []);
        continue
    end

    producedSubset = combinedEvents(idx);
    psth(sp.spike_times, [producedSubset.t_on], cfg.produced_window_s, cfg.dt, ...
        'Axes', ax, 'Title', producedSets{kk}.title);
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
