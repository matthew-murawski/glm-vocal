function [Xblk, info] = build_history_block(sps, stim, window_s)
% section defaults and input shaping
% we supply a default history window and reshape the spike counts to a column vector so the later sparse assembly stays simple.
mustHaveField(stim, 'dt');
if nargin < 3 || isempty(window_s)
    window_s = [stim.dt, stim.dt];
end
sps = double(sps(:));

% section input checks
% we verify the spike raster shape, ensure the stimulus struct has the required fields, and confirm the requested window is compatible with the timebase.
mustHaveField(stim, 't');
nT = numel(stim.t);
if numel(sps) ~= nT
    error('build_history_block:SizeMismatch', 'spike raster length (%d) must match stim.t (%d).', numel(sps), nT);
end
if ~isnumeric(window_s) || numel(window_s) ~= 2
    error('build_history_block:WindowShape', 'window_s must be a numeric two-element vector.');
end
window_s = double(window_s(:)');
if window_s(1) > window_s(2)
    error('build_history_block:WindowOrder', 'window_s must have a nondecreasing start and stop.');
end

dt = stim.dt;
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('build_history_block:InvalidDt', 'stim.dt must be a positive scalar.');
end

lagStartBins = round(window_s(1) / dt);
lagEndBins = round(window_s(2) / dt);
lagTol = 1e-9;
if abs(window_s(1) - lagStartBins * dt) > lagTol || abs(window_s(2) - lagEndBins * dt) > lagTol
    error('build_history_block:WindowResolution', 'window_s bounds must align with integer multiples of stim.dt.');
end
if lagStartBins <= 0
    error('build_history_block:HistoryWindow', 'history window must exclude lag 0; set start ≥ stim.dt.');
end

lagOffsets = lagStartBins:lagEndBins;
nLags = numel(lagOffsets);

% section sparse assembly
% we gather the valid history samples for each lag, skipping bins without sufficient lookback so they stay implicitly zero.
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
    valIdx{c} = sps(src(valid));
end

rows = vertcat(rowIdx{:});
cols = vertcat(colIdx{:});
vals = vertcat(valIdx{:});

if isempty(rows)
    Xblk = sparse(nT, nLags);
else
    Xblk = sparse(rows, cols, vals, nT, nLags);
end

% section metadata
% we expose the lag bins and times so downstream packing logic can keep track of the spike-history parameters.
info = struct();
info.mode = 'history';
info.window_s = window_s;
info.lag_bins = lagOffsets(:);
info.lag_times_s = info.lag_bins * dt;
end

function mustHaveField(s, fname)
% section struct guard
% small helper to ensure required stimulus fields are present before we proceed with computation.
if ~isstruct(s) || ~isfield(s, fname)
    error('build_history_block:MissingField', 'stim.%s is required.', fname);
end
end
