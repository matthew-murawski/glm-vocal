function fig = plot_predictions(stim, sps, rate, ev, start_sec, duration_sec, plot_title)
% section input checks
% confirm basic alignment and parameters before attempting to render the comparison plot.
if nargin < 7
    error('plot_predictions:NotEnoughInputs', 'plot_predictions requires seven input arguments.');
end
if ~isstruct(stim) || ~isfield(stim, 't') || ~isfield(stim, 'dt')
    error('plot_predictions:BadStim', 'stim must be a struct with fields t and dt.');
end
if ~isnumeric(stim.dt) || ~isscalar(stim.dt) || ~isfinite(stim.dt) || stim.dt <= 0
    error('plot_predictions:BadDt', 'stim.dt must be a positive scalar.');
end
t = stim.t(:);
sps = sps(:);
rate = rate(:);
if numel(t) ~= numel(sps) || numel(t) ~= numel(rate)
    error('plot_predictions:SizeMismatch', 'stim.t, sps, and rate must have matching lengths.');
end
if ~isnumeric(start_sec) || ~isscalar(start_sec) || ~isfinite(start_sec)
    error('plot_predictions:BadStart', 'start_sec must be a finite scalar.');
end
if ~isnumeric(duration_sec) || ~isscalar(duration_sec) || ~isfinite(duration_sec) || duration_sec <= 0
    error('plot_predictions:BadDuration', 'duration_sec must be a positive finite scalar.');
end
if ~ischar(plot_title) && ~isstring(plot_title)
    error('plot_predictions:BadTitle', 'plot_title must be a character vector or string scalar.');
end
plot_title = char(plot_title);
% section window selection
% choose the indices that fall inside the requested window, clipping to available samples where needed.
window_start = start_sec;
window_end = start_sec + duration_sec;
mask = (t >= window_start) & (t <= window_end);
if ~any(mask)
    error('plot_predictions:EmptyWindow', 'the requested window contains no samples.');
end
% section convert to rates
% compute firing-rate estimates (hz) so actual and predicted traces share units.
dt = stim.dt;
t_window = t(mask);
actual_hz = sps(mask) ./ dt;
pred_hz = rate(mask) ./ dt;
% section prepare events
% gather onset/offset pairs and kinds so call annotations become span boxes below the rate plot.
event_on = [];
event_off = [];
event_kinds = {};
if nargin >= 4 && ~isempty(ev) && isstruct(ev)
    has_on = all(isfield(ev, 't_on'));
    has_kind = all(isfield(ev, 'kind'));
    if has_on && has_kind
        for ii = 1:numel(ev)
            evt = ev(ii);
            on_val = evt.t_on;
            if ~isfinite(on_val)
                continue
            end
            if isfield(evt, 't_off')
                off_val = evt.t_off;
            else
                off_val = on_val;
            end
            if ~isfinite(off_val) || off_val < on_val
                off_val = on_val;
            end
            kind_val = evt.kind;
            if isstring(kind_val)
                kind_val = char(kind_val);
            end
            kind_norm = lower(strtrim(kind_val));
            if isempty(kind_norm)
                continue
            end
            if ~strcmp(kind_norm, 'produced') && ~strcmp(kind_norm, 'perceived')
                continue
            end
            event_on(end+1, 1) = on_val; %#ok<AGROW>
            event_off(end+1, 1) = off_val; %#ok<AGROW>
            event_kinds{end+1, 1} = kind_norm; %#ok<AGROW>
        end
    end
end
% section draw lines
% render the actual and predicted rate traces and seed legend entries for call markers.
actual_color = [0.1, 0.1, 0.1];
pred_color = [0.0, 0.6, 0.2];
produced_color = [0.85, 0.33, 0.1];
perceived_color = [0.0, 0.45, 0.74];
fig = figure('Color', 'w');
layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
ax_rate = nexttile(layout, [3, 1]);
hold(ax_rate, 'on');

% --- MODIFICATION IS HERE ---
actual_handle = stem(ax_rate, t_window, actual_hz, ...
    'LineStyle', '--', ...
    'Marker', '.', ...
    'MarkerSize', 12, ...
    'Color', actual_color, ...
    'DisplayName', 'Actual Rate', ...
    'BaseValue', 0);
% --- END MODIFICATION ---

pred_handle = plot(ax_rate, t_window, pred_hz, 'Color', pred_color, 'LineWidth', 1.5, 'DisplayName', 'Predicted Rate');
legend_handles = [actual_handle, pred_handle];
legend_labels = {'Actual Rate', 'Predicted Rate'};

overlaps_window = false(numel(event_on), 1);
for ii = 1:numel(event_on)
    overlaps_window(ii) = (event_off(ii) >= window_start) && (event_on(ii) <= window_end);
end

has_produced = any(overlaps_window & strcmp(event_kinds, 'produced'));
has_perceived = any(overlaps_window & strcmp(event_kinds, 'perceived'));
if has_produced
    legend_handles(end+1) = plot(ax_rate, nan, nan, 'Color', produced_color, 'LineWidth', 1.2, 'DisplayName', 'Produced Call'); %#ok<AGROW>
    legend_labels{end+1} = 'Produced Call'; %#ok<AGROW>
end
if has_perceived
    legend_handles(end+1) = plot(ax_rate, nan, nan, 'Color', perceived_color, 'LineWidth', 1.2, 'DisplayName', 'Perceived Call'); %#ok<AGROW>
    legend_labels{end+1} = 'Perceived Call'; %#ok<AGROW>
end

ax_events = nexttile(layout);
hold(ax_events, 'on');
for ii = 1:numel(event_on)
    if ~overlaps_window(ii)
        continue
    end
    span_start = max(event_on(ii), window_start);
    span_end = min(event_off(ii), window_end);
    if span_end <= span_start
        continue
    end
    kind = event_kinds{ii};
    if strcmp(kind, 'produced')
        y_bottom = 0.55;
        y_top = 0.95;
        color_val = produced_color;
    else
        y_bottom = 0.05;
        y_top = 0.45;
        color_val = perceived_color;
    end
    patch(ax_events, ...
        'XData', [span_start, span_start, span_end, span_end], ...
        'YData', [y_bottom, y_top, y_top, y_bottom], ...
        'FaceColor', 'none', ...
        'EdgeColor', color_val, ...
        'LineWidth', 1.4);
end
hold(ax_events, 'off');

hold(ax_rate, 'off');
% section finalize plot
% polish axis limits, labels, grid, and title so the plot communicates context immediately.
xlim(ax_rate, [window_start, window_end]);
xlim(ax_events, [window_start, window_end]);
ylabel(ax_rate, 'Firing Rate (Hz)');
title(ax_rate, plot_title);
legend(ax_rate, legend_handles, legend_labels, 'Location', 'best');
grid(ax_rate, 'on');

ylim(ax_events, [0, 1]);
yticks(ax_events, [0.25, 0.75]);
yticklabels(ax_events, {'perceived', 'produced'});
xlabel(ax_events, 'Time (s)');
ylabel(ax_events, 'Events');
grid(ax_events, 'off');
end
