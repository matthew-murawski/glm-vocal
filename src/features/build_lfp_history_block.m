function [Xblk, info] = build_lfp_history_block(lfp_signal, stim, window_s, use_basis)
%BUILD_LFP_HISTORY_BLOCK Build autoregressive LFP history features.
%   [Xblk, info] = BUILD_LFP_HISTORY_BLOCK(lfp_signal, stim, window_s, use_basis)
%   creates lagged copies of the LFP signal for autoregressive modeling.
%
%   Inputs:
%     lfp_signal - (n_bins × 1) continuous LFP values aligned to stim.t
%     stim       - struct with fields .t and .dt
%     window_s   - [start, end] in seconds (e.g., [0.01, 0.5])
%                  must be causal (start > 0) to exclude lag 0
%     use_basis  - (optional) if true, project lags onto raised cosine basis
%
%   Returns:
%     Xblk - sparse (n_bins × n_lags) or (n_bins × n_basis) matrix
%     info - metadata struct with lag times and mode

% section defaults and input shaping
% we supply a default history window and reshape the lfp signal to a column vector
mustHaveField(stim, 'dt');
if nargin < 3 || isempty(window_s)
    window_s = [stim.dt, stim.dt];
end
if nargin < 4 || isempty(use_basis)
    use_basis = false;
end
lfp_signal = double(lfp_signal(:));

% section input checks
% we verify the LFP signal shape, ensure the stimulus struct has the required fields,
% and confirm the requested window is compatible with the timebase.
mustHaveField(stim, 't');
nT = numel(stim.t);
if numel(lfp_signal) ~= nT
    error('build_lfp_history_block:SizeMismatch', ...
          'LFP signal length (%d) must match stim.t (%d).', numel(lfp_signal), nT);
end
if ~isnumeric(window_s) || numel(window_s) ~= 2
    error('build_lfp_history_block:WindowShape', 'window_s must be a numeric two-element vector.');
end
window_s = double(window_s(:)');
if window_s(1) > window_s(2)
    error('build_lfp_history_block:WindowOrder', 'window_s must have a nondecreasing start and stop.');
end

dt = stim.dt;
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('build_lfp_history_block:InvalidDt', 'stim.dt must be a positive scalar.');
end

lagStartBins = round(window_s(1) / dt);
lagEndBins = round(window_s(2) / dt);
lagTol = 1e-9;
if abs(window_s(1) - lagStartBins * dt) > lagTol || abs(window_s(2) - lagEndBins * dt) > lagTol
    error('build_lfp_history_block:WindowResolution', ...
          'window_s bounds must align with integer multiples of stim.dt.');
end
if lagStartBins <= 0
    error('build_lfp_history_block:HistoryWindow', ...
          'LFP history window must exclude lag 0; set start ≥ stim.dt.');
end

lagOffsets = lagStartBins:lagEndBins;
nLags = numel(lagOffsets);

% section sparse assembly
% we gather the valid LFP history samples for each lag, skipping bins without
% sufficient lookback so they stay implicitly zero.
rowIdx = cell(nLags, 1);
colIdx = cell(nLags, 1);
valIdx = cell(nLags, 1);

allRows = (1:nT)';
for c = 1:nLags
    lag = lagOffsets(c);
    src = allRows - lag;
    valid = src >= 1 & src <= nT;

    rowIdx{c} = allRows(valid);
    colIdx{c} = repmat(c, nnz(valid), 1);
    valIdx{c} = lfp_signal(src(valid));
end

rows = vertcat(rowIdx{:});
cols = vertcat(colIdx{:});
vals = vertcat(valIdx{:});

if isempty(rows)
    Xblk_raw = sparse(nT, nLags);
else
    Xblk_raw = sparse(rows, cols, vals, nT, nLags);
end

% section optional basis projection
% if requested, project the raw lagged LFP onto a smooth raised cosine basis
% to reduce dimensionality and enforce smoothness
if use_basis && nLags > 1
    % build raised cosine basis with default parameters
    n_basis = min(8, nLags);  % use up to 8 basis functions
    overlap = 2.0;

    % create basis matrix (n_lags × n_basis)
    lag_times = lagOffsets * dt;
    B = build_raised_cosine_basis(lag_times, n_basis, overlap);

    % project: Xblk = Xblk_raw * B
    Xblk = Xblk_raw * B;

    info_mode = 'lfp_history_basis';
    info_basis = B;
else
    Xblk = Xblk_raw;
    info_mode = 'lfp_history';
    info_basis = [];
end

% section metadata
% we expose the lag bins and times so downstream packing logic can keep track
% of the LFP history parameters.
info = struct();
info.mode = info_mode;
info.window_s = window_s;
info.lag_bins = lagOffsets(:);
info.lag_times_s = lagOffsets(:) * dt;
if ~isempty(info_basis)
    info.basis = info_basis;
end

end

function mustHaveField(s, fname)
% section struct guard
% small helper to ensure required stimulus fields are present before we proceed with computation.
if ~isstruct(s) || ~isfield(s, fname)
    error('build_lfp_history_block:MissingField', 'stim.%s is required.', fname);
end
end

function B = build_raised_cosine_basis(lag_times, n_basis, overlap)
% BUILD_RAISED_COSINE_BASIS Create a raised cosine basis for smoothing.
%   B = BUILD_RAISED_COSINE_BASIS(lag_times, n_basis, overlap)
%   creates a matrix of raised cosine functions.
%
%   Inputs:
%     lag_times - vector of lag times in seconds
%     n_basis   - number of basis functions
%     overlap   - overlap factor (controls width of each cosine)
%
%   Returns:
%     B - (n_lags × n_basis) basis matrix

n_lags = numel(lag_times);
t_min = min(lag_times);
t_max = max(lag_times);
t_range = t_max - t_min;

if t_range == 0
    % degenerate case: single lag
    B = ones(n_lags, 1);
    return;
end

% create evenly spaced centers across the lag window
centers = linspace(t_min, t_max, n_basis);
span = overlap * t_range / (n_basis - 1);

% build each basis function as a raised cosine
B = zeros(n_lags, n_basis);
for k = 1:n_basis
    center = centers(k);
    for i = 1:n_lags
        t = lag_times(i);
        if abs(t - center) <= span
            B(i, k) = 0.5 * (cos(pi * (t - center) / span) + 1);
        end
    end

    % L1 normalization (sum to 1)
    col_sum = sum(B(:, k));
    if col_sum > 0
        B(:, k) = B(:, k) / col_sum;
    end
end

end
