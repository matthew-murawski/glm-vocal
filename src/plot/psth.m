function out = psth(spike_times, event_times, win, bin_size, varargin)
% make a psth figure using align_and_bin_spikes.
% supports 1..N event-time sets (cell array)
%
% usage:
% out = psth(spk, ev, [-1 1], 0.02);
% out = psth(spk, {evA, evB, evC}, [-1 1], 0.02, 'Labels',{'heard','silence','produced'});
% out = psth(spk, ev, [-1 1], 0.02, 'SmoothingBins', 3, 'PlotSEM', true);
% out = psth(spk, ev, [-1 1], 0.02, 'Title','unit A4 twitter trials');
% out = psth(spk, ev, [-1 1], 0.02, 'XLabel','time (s) relative to call onset');
% out = psth(spk, ev, [-1 1], 0.02, 'ZScore', true); % per-trial z across bins
%
% inputs
%   spike_times : vector of spike times (s)
%   event_times : vector of event times (s) or a cell {ev1, ev2, ..., evN}
%   win         : [t_pre t_post] window around each event (s), e.g., [-1 1]
%   bin_size    : bin width (s), e.g., 0.02 for 20 ms
%
% name–value options (all optional)
%   'Axes'          : existing axes handle to plot into (default: new figure)
%   'SmoothingBins' : movmean window (bins) applied to mean psth (default: 3; set 1 for none)
%   'PlotSEM'       : true/false to draw ±SEM shaded region (default: true)
%   'LineWidth'     : line width for mean trace (default: 2)
%   'ShadeAlpha'    : facealpha for sem patch (default: 0.3)
%   'Color'         : rgb triple or K×3 matrix for K traces (default: colororder)
%   'Labels'        : labels for legend (cellstr or string), length ≤ K (default: {})
%   'Title'         : custom title string; empty → auto title (default: '')
%   'XLabel'        : custom x-axis label (default: 'time from event (s)')
%   'ZScore'        : if true, z-scores counts per trial across bins (default: false)
%
% outputs (struct)
%   out.conditions(k): per-condition results: counts, rates, psth_mean, psth_sem
%   out.bin_edges     : shared edges (1×(bins+1))
%   out.bin_centers   : shared centers (1×bins)
%   out.fig           : figure handle (if created here)
%   out.ax            : axes handle used
%   out.h_mean        : line handles (1×K)
%   out.h_sem         : patch handles (1×K; [] if not drawn)
%   out.is_zscored    : true/false
%   out.y_units       : 'spikes/s' or 'z-score'
%
%   for backward-compatibility, condition 1 is mirrored to:
%     out.counts, out.rates, out.psth_mean, out.psth_sem

% parse options
ip = inputParser;
ip.addParameter('Axes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
ip.addParameter('SmoothingBins', 3, @(x) isscalar(x) && x >= 1);
ip.addParameter('PlotSEM', true, @(x) islogical(x) || isnumeric(x));
ip.addParameter('LineWidth', 2, @(x) isscalar(x) && x > 0);
ip.addParameter('ShadeAlpha', 0.3, @(x) isscalar(x) && x >= 0 && x <= 1);
ip.addParameter('Color', [], @(x) isempty(x) || (isnumeric(x) && (numel(x)==3 || size(x,2)==3)));
ip.addParameter('Labels', {}, @(x) iscellstr(x) || isstring(x) || isempty(x));
ip.addParameter('Title', '', @(x) ischar(x) || isstring(x));
ip.addParameter('XLabel', 'time from event (s)', @(x) ischar(x) || isstring(x));
ip.addParameter('ZScore', false, @(x) islogical(x) || isnumeric(x));
ip.parse(varargin{:});
opt = ip.Results;
opt.ZScore = logical(opt.ZScore);
labels_in = cellstr(opt.Labels);

% normalize event_times to a cell for uniform logic
if iscell(event_times)
    ets = event_times(:);
else
    ets = {event_times};
end
K = numel(ets);
if K < 1
    error('event_times must be a vector or a nonempty cell array of vectors.');
end

% compute per-condition aligned counts and summary stats
conds = repmat(struct( ...
    'counts', [], 'rates', [], 'psth_mean', [], 'psth_sem', [], ...
    'n_trials', 0, 'bin_edges', [], 'bin_centers', []), 1, K);

for k = 1:K
    evk = ets{k};
    [counts_k, edges_k, centers_k] = align_and_bin_spikes( ...
        spike_times, evk, win, bin_size, 'ZScore', opt.ZScore);

    n_trials_k = size(counts_k, 1);
    n_bins_k   = size(counts_k, 2);

    if n_trials_k == 0 || n_bins_k == 0
        rates_k     = [];
        psth_mean_k = [];
        psth_sem_k  = [];
    else
        if opt.ZScore
            rates_k = counts_k;            % already unitless
        else
            rates_k = counts_k / bin_size; % spikes/s
        end
        psth_mean_k = mean(rates_k, 1);
        if n_trials_k > 1
            psth_sem_k = std(rates_k, 0, 1) / sqrt(n_trials_k);
        else
            psth_sem_k = zeros(1, n_bins_k);
        end
    end

    conds(k).counts       = counts_k;
    conds(k).rates        = rates_k;
    conds(k).psth_mean    = psth_mean_k;
    conds(k).psth_sem     = psth_sem_k;
    conds(k).n_trials     = n_trials_k;
    conds(k).bin_edges    = edges_k;
    conds(k).bin_centers  = centers_k;
end

% pick units for labeling
y_units = tern(opt.ZScore, 'z-score', 'spikes/s');

% ensure identical bins across conditions
ref_edges = conds(1).bin_edges;
for k = 2:K
    be = conds(k).bin_edges;
    if ~isequal(size(ref_edges), size(be)) || any(abs(ref_edges - be) > eps(max(1,abs(ref_edges))))
        error('bin edges differ between conditions; ensure identical win/bin_size.');
    end
end

% early exit if nothing to plot
if all(arrayfun(@(c) isempty(c.rates) || isempty(c.psth_mean), conds))
    warning('no trials or no bins; nothing to plot.');
    out = struct( ...
        'conditions', conds, ...
        'counts', conds(1).counts, ...
        'rates',  conds(1).rates, ...
        'psth_mean', conds(1).psth_mean, ...
        'psth_sem',  conds(1).psth_sem, ...
        'bin_edges', conds(1).bin_edges, ...
        'bin_centers', conds(1).bin_centers, ...
        'fig', [], 'ax', [], 'h_mean', [], 'h_sem', [], ...
        'is_zscored', opt.ZScore, 'y_units', y_units);
    return
end

% set up axes/figure and colors
made_fig = false;
if isempty(opt.Axes)
    fig = figure('Color','w');
    ax  = axes(fig); %#ok<LAXES>
    made_fig = true;
else
    ax  = opt.Axes;
    fig = ancestor(ax, 'figure');
    hold(ax, 'on');
end
ax.LineWidth = 0.5;
hold(ax, 'on');

colors = choose_colors(ax, opt.Color, K);

% draw sem patches first, then mean traces
h_sem  = gobjects(1, K);
h_mean = gobjects(1, K);
for k = 1:K
    if isempty(conds(k).psth_mean)
        h_sem(k)  = gobjects(1);
        h_mean(k) = gobjects(1);
        continue
    end

    mu = conds(k).psth_mean;
    se = conds(k).psth_sem;

    % apply the same smoothing to both mean and sem for plotting
    if opt.SmoothingBins > 1 && ~isempty(mu)
        mu_plot = movmean(mu, opt.SmoothingBins);
        se_plot = movmean(se, opt.SmoothingBins);
    else
        mu_plot = mu;
        se_plot = se;
    end

    h_sem(k) = gobjects(1);
    if islogical(opt.PlotSEM) && opt.PlotSEM && any(se_plot > 0)
        x = [conds(k).bin_centers, fliplr(conds(k).bin_centers)];
        y = [mu_plot + se_plot,    fliplr(mu_plot - se_plot)];
        h_sem(k) = fill(ax, x, y, colors(k,:), 'FaceAlpha', opt.ShadeAlpha, 'EdgeColor','none');
    end

    h_mean(k) = plot(ax, conds(k).bin_centers, mu_plot, ...
        'Color', colors(k,:), 'LineWidth', opt.LineWidth);
end

% axis cosmetics
xline(ax, 0, ':', 'Color', [0 0 0 0.5], 'LineWidth', 1);
xlabel(ax, char(opt.XLabel));
ylabel(ax, y_units);
box(ax, 'off');
set(ax, 'Layer','top');

% legend
valid = arrayfun(@(h) isgraphics(h), h_mean);
idx   = find(valid);
if K > 1 && ~isempty(labels_in)
    m = min(numel(labels_in), numel(idx));
    if m >= 2
        legend(ax, h_mean(idx(1:m)), labels_in(1:m), 'Box','off', 'Location','best');
    end
end

% title
if strlength(string(opt.Title)) > 0
    ttl = char(opt.Title);
else
    nvec = arrayfun(@(c) c.n_trials, conds);
    if opt.ZScore
        ttl = sprintf('psth (z-scored, n per cond = [%s], bin=%.0f ms)', ...
            strjoin(string(nvec), ' '), bin_size*1e3);
    else
        ttl = sprintf('psth (n per cond = [%s], bin=%.0f ms)', ...
            strjoin(string(nvec), ' '), bin_size*1e3);
    end
end
title(ax, ttl);

% assemble outputs
out = struct();
out.conditions   = conds;
out.counts       = conds(1).counts;
out.rates        = conds(1).rates;
out.psth_mean    = conds(1).psth_mean;
out.psth_sem     = conds(1).psth_sem;
out.bin_edges    = conds(1).bin_edges;
out.bin_centers  = conds(1).bin_centers;
out.fig          = [];
out.ax           = ax;
out.h_mean       = h_mean;
out.h_sem        = h_sem;
out.is_zscored   = opt.ZScore;
out.y_units      = y_units;
if made_fig, out.fig = fig; end
end

% tiny helper to keep code tidy
function y = tern(cond, a, b)
if cond, y = a; else, y = b; end
end

% choose k colors, respecting any provided color(s)
function C = choose_colors(ax, color_opt, K)
co = get(ax, 'ColorOrder');
if isempty(color_opt)
    C = zeros(K,3);
    start = numel(ax.Children);  % rough cycling
    for k = 1:K
        idx = mod(start + k - 1, size(co,1)) + 1;
        C(k,:) = co(idx,:);
    end
else
    color_opt = double(color_opt);
    if isvector(color_opt) && numel(color_opt)==3
        C = repmat(color_opt(:).', K, 1);
    elseif size(color_opt,2) == 3
        C = color_opt(1:min(size(color_opt,1),K), :);
        if size(C,1) < K
            % pad by cycling colororder if fewer rows than K
            pad = zeros(K-size(C,1),3);
            for i = 1:size(pad,1)
                pad(i,:) = co(mod(i-1, size(co,1)) + 1, :);
            end
            C = [C; pad];
        end
    else
        error('Color must be 1×3 or K×3 (K = number of traces).');
    end
end
end