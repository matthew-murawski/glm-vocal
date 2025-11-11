function binned = bin_high_gamma(hg, stim, cfg)
%BIN_HIGH_GAMMA Downsample/aggregate high gamma power to GLM time base.
%   binned = BIN_HIGH_GAMMA(hg, stim, cfg) aggregates continuous high gamma
%   power trace to match the time bins defined in stim.t.
%
%   Inputs:
%       hg      - struct from load_high_gamma with fields:
%                 * power: [nSamples × 1] high gamma power trace
%                 * fs: sampling rate (Hz)
%                 * t: [nSamples × 1] time vector (seconds)
%       stim    - struct with field:
%                 * t: [nBins × 1] GLM time bin centers
%       cfg     - struct with field:
%                 * preprocess.hg_binning_method: 'mean'|'median'|'rms'|'max'
%
%   Output:
%       binned  - struct with fields:
%                 * power: [nBins × 1] aggregated power
%                 * t: stim.t (for reference)
%                 * method: aggregation method used
%                 * n_empty_bins: count of interpolated bins
%                 * warnings: cell array of warning messages

% validate inputs
if ~isstruct(hg) || ~all(isfield(hg, {'power', 'fs', 't'}))
    error('glm:InvalidInput', 'hg must be a struct with fields: power, fs, t');
end
if ~isstruct(stim) || ~isfield(stim, 't')
    error('glm:InvalidInput', 'stim must be a struct with field: t');
end
if ~isstruct(cfg) || ~isfield(cfg, 'preprocess') || ~isfield(cfg.preprocess, 'hg_binning_method')
    error('glm:InvalidInput', 'cfg must have cfg.preprocess.hg_binning_method');
end

% extract parameters
power_in = hg.power(:);
t_in = hg.t(:);
t_out = stim.t(:);
method = cfg.preprocess.hg_binning_method;
fs = hg.fs;

% validate method
valid_methods = {'mean', 'median', 'rms', 'max'};
if ~ismember(lower(method), valid_methods)
    error('glm:InvalidInput', 'hg_binning_method must be one of: %s', strjoin(valid_methods, ', '));
end

nBinsOut = length(t_out);
nSamplesIn = length(power_in);

% initialize output and tracking
power_out = nan(nBinsOut, 1);
empty_bins = false(nBinsOut, 1);
warnings_list = {};

% compute bin width from stim.t
if nBinsOut > 1
    dt = median(diff(t_out));
else
    dt = 0.01;  % default
end

% check for upsampling condition
if fs < (1 / dt)
    warnings_list{end+1} = sprintf(...
        'High gamma sampling rate (%.2f Hz) is lower than bin rate (%.2f Hz). Consider using dt >= %.4f s.', ...
        fs, 1/dt, 1/fs);
end

% bin boundaries: each bin centered at t_out(i) spans [t_out(i) - dt/2, t_out(i) + dt/2)
% for efficiency, we'll iterate through bins and find samples in each
for i = 1:nBinsOut
    % define bin edges
    if i == 1
        % first bin: left edge at t_out(1) - dt/2 or min(t_in), whichever is larger
        left_edge = max(t_out(i) - dt/2, min(t_in) - eps);
    else
        % use midpoint between current and previous bin
        left_edge = (t_out(i-1) + t_out(i)) / 2;
    end

    if i == nBinsOut
        % last bin: right edge at t_out(end) + dt/2 or max(t_in), whichever is smaller
        right_edge = min(t_out(i) + dt/2, max(t_in) + eps);
    else
        % use midpoint between current and next bin
        right_edge = (t_out(i) + t_out(i+1)) / 2;
    end

    % find all samples in this bin
    in_bin = (t_in >= left_edge) & (t_in < right_edge);
    samples_in_bin = power_in(in_bin);

    if isempty(samples_in_bin)
        % no samples in this bin - mark for interpolation
        empty_bins(i) = true;
    else
        % aggregate samples using specified method
        power_out(i) = aggregate_samples(samples_in_bin, method);
    end
end

% handle empty bins via interpolation
n_empty = sum(empty_bins);
if n_empty > 0
    warnings_list{end+1} = sprintf(...
        'Found %d empty bins (%.2f%%). Using linear interpolation.', ...
        n_empty, 100 * n_empty / nBinsOut);

    power_out = interpolate_empty_bins(power_out, empty_bins);
end

% ensure non-negative output
if any(power_out < 0)
    n_negative = sum(power_out < 0);
    warnings_list{end+1} = sprintf(...
        'Binning produced %d negative values. Setting to zero.', n_negative);
    power_out(power_out < 0) = 0;
end

% build output struct
binned = struct();
binned.power = power_out;
binned.t = t_out;
binned.method = method;
binned.n_empty_bins = n_empty;
binned.warnings = warnings_list;

end

%% helper functions

function result = aggregate_samples(samples, method)
%AGGREGATE_SAMPLES Apply aggregation method to samples within a bin.

switch lower(method)
    case 'mean'
        result = mean(samples);
    case 'median'
        result = median(samples);
    case 'rms'
        result = sqrt(mean(samples.^2));
    case 'max'
        result = max(samples);
    otherwise
        error('glm:InvalidMethod', 'Unknown aggregation method: %s', method);
end
end

function power_out = interpolate_empty_bins(power_in, empty_idx)
%INTERPOLATE_EMPTY_BINS Fill empty bins using linear interpolation.
%   Handles edge cases where leading or trailing bins are empty.

power_out = power_in;
n = length(power_in);

% find indices of non-empty bins
filled_idx = find(~empty_idx);

if isempty(filled_idx)
    % all bins are empty - cannot interpolate, set to zero
    power_out(:) = 0;
    return
end

% handle leading empty bins (before first filled bin)
if filled_idx(1) > 1
    % extrapolate using first valid value
    power_out(1:filled_idx(1)-1) = power_in(filled_idx(1));
end

% handle trailing empty bins (after last filled bin)
if filled_idx(end) < n
    % extrapolate using last valid value
    power_out(filled_idx(end)+1:n) = power_in(filled_idx(end));
end

% handle interior empty bins using linear interpolation
empty_interior = find(empty_idx);
for i = 1:length(empty_interior)
    idx = empty_interior(i);

    % skip if leading or trailing (already handled)
    if idx < filled_idx(1) || idx > filled_idx(end)
        continue
    end

    % find nearest filled bins on left and right
    left_idx = filled_idx(find(filled_idx < idx, 1, 'last'));
    right_idx = filled_idx(find(filled_idx > idx, 1, 'first'));

    if ~isempty(left_idx) && ~isempty(right_idx)
        % linear interpolation
        alpha = (idx - left_idx) / (right_idx - left_idx);
        power_out(idx) = (1 - alpha) * power_in(left_idx) + alpha * power_in(right_idx);
    end
end

end
