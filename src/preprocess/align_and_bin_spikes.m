function [counts, bin_edges, bin_centers] = align_and_bin_spikes(spike_times, event_times, win, bin_size, varargin)
% align spikes to events and bin them.

% inputs
% spike_times: vector of spike times in seconds.
% event_times: vector of event times in seconds.
% win: 1x2 time window around each event, e.g. [-0.1, 1]
% bin_size: bin width in seconds, e.g. 0.02 for 20ms bins
% name-value:
%   'ZScore' (false)  -> if true, z-scores counts per trial across bins

% outputs
% counts: [n_events x n_bins] spike counts per trial and bin (z-scored if requested).
% bin_edges: 1 x (n_bins+1) vector of bin edges (relative time).
% bin_centers: 1 x n_bins vector of bin centers (relative time).

    %% parse optional args
    zscore_flag = false;
    if ~isempty(varargin)
        if mod(numel(varargin),2) ~= 0
            error('name-value pairs must come in pairs.');
        end
        for i = 1:2:numel(varargin)
            name = varargin{i};
            val  = varargin{i+1};
            if ~(ischar(name) || isstring(name))
                error('parameter names must be char or string.');
            end
            switch lower(char(name))
                case 'zscore'
                    zscore_flag = logical(val);
                otherwise
                    error('unrecognized parameter name: %s', char(name));
            end
        end
    end

    %% basic
    if isempty(spike_times)
        spike_times = [];
    end
    if isempty(event_times)
        counts = zeros(0,0);
        bin_edges = [];
        bin_centers = [];
        return
    end

    % make em column vectors
    spike_times = spike_times(:);
    event_times = event_times(:);

    %% construct bin edges
    % regular edges with a hard stop at t1 (right edge)
    t0 = win(1);
    t1 = win(2);
    if ~(isscalar(bin_size) && bin_size > 0)
        error('bin_size must be a positive scalar (seconds).');
    end
    if t1 <= t0
        error('win must be [t_pre t_post] with t_post > t_pre.');
    end
    bin_edges = t0:bin_size:t1;
    if bin_edges(end) < t1
        bin_edges = [bin_edges, t1];
    else
        bin_edges(end) = t1; % clamp tiny float drift
    end
    n_bins = numel(bin_edges) - 1;
    n_ev   = numel(event_times);
    counts = zeros(n_ev, n_bins);

    %% core loop: align spikes to each event and bin
    % left-inclusive, right-exclusive bins match histcounts defaults
    for k = 1:n_ev
        rel = spike_times - event_times(k);
        in  = (rel >= t0) & (rel < t1);
        if any(in)
            counts(k, :) = histcounts(rel(in), bin_edges);
        else
            counts(k, :) = 0;
        end
    end

    %% optional per-trial zscoring across bins
    % this removes trial-wise baseline offsets without needing a separate baseline window
    if zscore_flag && n_ev > 0 && n_bins > 0
        mu = mean(counts, 2);
        sd = std(counts, 0, 2);
        sd_safe = sd;
        sd_safe(sd_safe == 0) = 1; % avoid divide-by-zero
        counts = bsxfun(@minus, counts, mu);
        counts = bsxfun(@rdivide, counts, sd_safe);
        zero_var_rows = (sd == 0);
        if any(zero_var_rows)
            counts(zero_var_rows, :) = 0; % define z-score as 0 when variance is 0
        end
    end

    %% bin centers for plotting
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
end
