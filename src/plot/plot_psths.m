function plot_psths(sp, heardEvents, producedEvents, cfg, outdir)
% makes a figure with two psth subplots: one for heard, one for produced.
% this helps compare the raw firing rate patterns to the glm kernels.

% section 1: set up the figure
% we'll create a figure to hold the two subplots. this is created with
% 'Visible','off' because this is intended for a script, not interactive use.
fig = figure('Name', 'PSTHs', 'Color', 'w', 'Position', [100, 100, 800, 600], 'Visible', 'off');

% section 2: plot psth for heard calls
% this subplot shows the neuron's response aligned to the onset of calls it heard.
ax1 = subplot(2, 1, 1);
if ~isempty(heardEvents)
    % we use the existing psth.m function to do the heavy lifting.
    % the window and bin size are taken from the config to match the glm.
    psth(sp.spike_times, [heardEvents.t_on], cfg.heard_window_s, cfg.dt, ...
        'Axes', ax1, 'Title', 'psth aligned to heard calls');
else
    % if there are no heard events, we just show an empty plot with a title.
    title(ax1, 'psth aligned to heard calls (no events)');
    box(ax1, 'on');
    set(ax1, 'xtick', [], 'ytick', []);
end

% section 3: plot psth for produced calls
% this subplot shows the neuron's response aligned to its own calls.
ax2 = subplot(2, 1, 2);
if ~isempty(producedEvents)
    % same as above, but for the produced calls.
    % note that the window for produced calls can be different from heard calls.
    psth(sp.spike_times, [producedEvents.t_on], cfg.produced_window_s, cfg.dt, ...
        'Axes', ax2, 'Title', 'psth aligned to produced calls');
else
    % handle the case with no produced events.
    title(ax2, 'psth aligned to produced calls (no events)');
    box(ax2, 'on');
    set(ax2, 'xtick', [], 'ytick', []);
end

% section 4: save the figure
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