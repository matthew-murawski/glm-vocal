function [Xblk, info] = build_kernel_block(z, stim, window_s, mode)
% section defaults and input shaping
% we supply defaults for optional arguments and reshape the regressor so the rest of the logic can rely on column-oriented data. this keeps downstream operations predictable.
if nargin < 3 || isempty(window_s)
    window_s = [0, 0];
end
if nargin < 4 || isempty(mode)
    mode = 'causal';
end
z = double(z(:));

% section input checks
% we verify that the regressor matches the stimulus grid and that the requested window can be represented on that grid. this keeps the sparse assembly free from hidden edge cases.
mustHaveField(stim, 'dt');
mustHaveField(stim, 't');
nT = numel(stim.t);
if numel(z) ~= nT
    error('build_kernel_block:SizeMismatch', 'regressor length (%d) must match stim.t (%d).', numel(z), nT);
end
if ~isnumeric(window_s) || numel(window_s) ~= 2
    error('build_kernel_block:WindowShape', 'window_s must be a numeric two-element vector.');
end
window_s = double(window_s(:)');
if window_s(1) > window_s(2)
    error('build_kernel_block:WindowOrder', 'window_s must have a nondecreasing start and stop.');
end
if ~ismember(mode, {'causal', 'symmetric'})
    error('build_kernel_block:ModeUnsupported', 'mode ''%s'' is not supported.', mode);
end

dt = stim.dt;
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('build_kernel_block:InvalidDt', 'stim.dt must be a positive scalar.');
end

lagStartBins = round(window_s(1) / dt);
lagEndBins = round(window_s(2) / dt);
lagTol = 1e-9;
if abs(window_s(1) - lagStartBins * dt) > lagTol || abs(window_s(2) - lagEndBins * dt) > lagTol
    error('build_kernel_block:WindowResolution', 'window_s bounds must align with integer multiples of stim.dt.');
end
if strcmp(mode, 'causal') && lagStartBins < 0
    error('build_kernel_block:CausalWindow', 'causal mode requires nonnegative window bounds.');
end
if strcmp(mode, 'symmetric')
    if lagStartBins > 0 || lagEndBins < 0
        error('build_kernel_block:SymmetricWindow', 'symmetric mode requires window_s to span negative through positive lags.');
    end
end

lagOffsets = lagStartBins:lagEndBins;
nLags = numel(lagOffsets);
if strcmp(mode, 'symmetric') && ~any(lagOffsets == 0)
    error('build_kernel_block:MissingZeroLag', 'symmetric mode requires window_s to include lag 0.');
end

% section sparse assembly
% we walk each lag, gather the overlapping samples, and populate the sparse design block with zero padding where the shift exceeds the available history. this keeps the block toeplitz-like without dense copies.
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
    valIdx{c} = z(src(valid));
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
% we return lag information so downstream code knows the temporal meaning of each column in the block. this makes later plotting and unpacking straightforward.
info = struct();
info.mode = mode;
info.window_s = window_s;
info.lag_bins = lagOffsets(:);
info.lag_times_s = info.lag_bins * dt;
end

function mustHaveField(s, fname)
% section struct guard
% this helper keeps the main function readable by checking that required stimulus fields exist. it fails fast when the stimulus struct is incomplete.
if ~isstruct(s) || ~isfield(s, fname)
    error('build_kernel_block:MissingField', 'stim.%s is required.', fname);
end
end
