function plot_predictions(stim, sps, rate, ev, start_sec, duration_sec, plot_title)
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
% gather event onset times and kinds so call annotations can be layered cleanly.
event_times = [];
event_kinds = {};
if nargin >= 4 && ~isempty(ev) && isstruct(ev) && isfield(ev, 't_on') && isfield(ev, 'kind')
    if numel(ev) > 1
        raw_times = [ev(:).t_on];
        raw_kinds = {ev(:).kind};
    else
        raw_times = ev.t_on;
        raw_kinds = ev.kind;
    end

    raw_times = raw_times(:);
    if iscell(raw_kinds)
        raw_kinds = raw_kinds(:);
    elseif isstring(raw_kinds)
        raw_kinds = cellstr(raw_kinds(:));
    elseif ischar(raw_kinds)
        if size(raw_kinds, 1) == numel(raw_times)
            raw_kinds = cellstr(raw_kinds);
        else
            raw_kinds = {raw_kinds};
        end
    else
        raw_kinds = cell(numel(raw_times), 1);
    end

    valid_times = [];
    valid_kinds = {};
    for ii = 1:numel(raw_times)
        time_val = raw_times(ii);
        if ~isfinite(time_val)
            continue
        end
        kind_val = '';
        if ii <= numel(raw_kinds) && ~isempty(raw_kinds{ii})
            kind_val = raw_kinds{ii};
        end
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
        valid_times(end+1, 1) = time_val; %#ok<AGROW>
        valid_kinds{end+1, 1} = kind_norm; %#ok<AGROW>
    end
    event_times = valid_times;
    event_kinds = valid_kinds;
end

% section draw lines
% render the actual and predicted rate traces and seed legend entries for call markers.
actual_color = [0.1, 0.1, 0.1];
pred_color = [0.0, 0.6, 0.2];
produced_color = [0.85, 0.33, 0.1];
perceived_color = [0.0, 0.45, 0.74];

fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
actual_handle = plot(ax, t_window, actual_hz, 'Color', actual_color, 'LineWidth', 1.5, 'DisplayName', 'Actual Rate');
pred_handle = plot(ax, t_window, pred_hz, 'Color', pred_color, 'LineWidth', 1.5, 'DisplayName', 'Predicted Rate');
legend_handles = [actual_handle, pred_handle];
legend_labels = {'Actual Rate', 'Predicted Rate'};

in_window_mask = (event_times >= window_start) & (event_times <= window_end);
has_produced = any(in_window_mask & strcmp(event_kinds, 'produced'));
has_perceived = any(in_window_mask & strcmp(event_kinds, 'perceived'));
if has_produced
    legend_handles(end+1) = plot(ax, nan, nan, 'Color', produced_color, 'LineWidth', 1.2, 'DisplayName', 'Produced Call'); %#ok<AGROW>
    legend_labels{end+1} = 'Produced Call'; %#ok<AGROW>
end
if has_perceived
    legend_handles(end+1) = plot(ax, nan, nan, 'Color', perceived_color, 'LineWidth', 1.2, 'DisplayName', 'Perceived Call'); %#ok<AGROW>
    legend_labels{end+1} = 'Perceived Call'; %#ok<AGROW>
end

% section event markers
% overlay vertical lines for produced and perceived onsets while keeping the legend tidy.
for ii = 1:numel(event_times)
    t_on = event_times(ii);
    if t_on < window_start || t_on > window_end
        continue
    end
    kind = event_kinds{ii};
    if strcmp(kind, 'produced')
        h = xline(ax, t_on, '-', 'Color', produced_color, 'LineWidth', 1.2);
        set(h, 'HandleVisibility', 'off');
    elseif strcmp(kind, 'perceived')
        h = xline(ax, t_on, '-', 'Color', perceived_color, 'LineWidth', 1.2);
        set(h, 'HandleVisibility', 'off');
    end
end
hold(ax, 'off');

% section finalize plot
% polish axis limits, labels, grid, and title so the plot communicates context immediately.
xlim(ax, [window_start, window_end]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Firing Rate (Hz)');
title(ax, plot_title);
legend(ax, legend_handles, legend_labels, 'Location', 'best');
grid(ax, 'on');
end
