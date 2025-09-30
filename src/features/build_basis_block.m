function [Xblk, info] = build_basis_block(z, stim, window_s, basisCfg)
% build a basis-projected design block from an impulse stream.
%
% z must already be aligned to the stimulus grid (column vector). window_s
% defines the symmetric lag window around the event, and basisCfg controls
% the smooth basis functions (currently raised cosines).

% section defaults and input checks
% we normalise optional inputs and reuse the kernel block helper to obtain lag-aligned copies.
if nargin < 4 || isempty(basisCfg)
    basisCfg = struct('kind', 'raised_cosine');
end

[lagBlk, lagInfo] = build_kernel_block(z, stim, window_s, 'symmetric');

kind = 'raised_cosine';
if isstruct(basisCfg) && isfield(basisCfg, 'kind') && ~isempty(basisCfg.kind)
    kind = char(basisCfg.kind);
end

switch kind
    case 'raised_cosine'
        [basisMatrix, basisMeta] = build_raised_cosine_basis(lagInfo.lag_times_s, basisCfg);
    otherwise
        error('build_basis_block:UnsupportedBasis', 'basis kind ''%s'' is not supported.', kind);
end

% section sparse assembly
% we project the lag block onto the basis functions to obtain smooth predictors.
Xdense = lagBlk * basisMatrix;
if issparse(Xdense)
    Xblk = Xdense;
else
    Xblk = sparse(Xdense);
end

% section metadata
% expose lag and basis details so downstream code can reconstruct kernels for plotting.
info = struct();
info.mode = basisMeta.mode;
info.window_s = window_s;
info.lag_bins = lagInfo.lag_bins;
info.lag_times_s = lagInfo.lag_times_s;
info.basis = basisMeta;
end

function [basisMatrix, meta] = build_raised_cosine_basis(lag_times_s, cfg)
% construct overlapping raised-cosine basis functions over the provided lag grid.

% section defaults
% we coerce configuration fields and derive basic spacing metrics.
if nargin < 2 || ~isstruct(cfg)
    cfg = struct();
end
if ~isfield(cfg, 'n_basis') || isempty(cfg.n_basis)
    cfg.n_basis = 6;
end
if ~isfield(cfg, 'overlap') || isempty(cfg.overlap)
    cfg.overlap = 2.0;
end
if ~isfield(cfg, 'normalize') || isempty(cfg.normalize)
    cfg.normalize = 'l1';
end

lag_times = double(lag_times_s(:));
nLags = numel(lag_times);
if nLags == 0
    error('build_raised_cosine_basis:EmptyGrid', 'lag_times_s must contain at least one entry.');
end

nBasis = double(cfg.n_basis);
if ~isscalar(nBasis) || ~isfinite(nBasis) || nBasis < 1 || floor(nBasis) ~= nBasis
    error('build_raised_cosine_basis:InvalidCount', 'n_basis must be a positive integer scalar.');
end
nBasis = int32(nBasis);

overlap = double(cfg.overlap);
if ~isscalar(overlap) || ~isfinite(overlap) || overlap <= 0
    error('build_raised_cosine_basis:InvalidOverlap', 'overlap must be a positive scalar.');
end

lagStep = 0;
if nLags > 1
    lagStep = median(diff(lag_times));
end
if lagStep <= 0
    lagStep = 1.0;
end

tMin = lag_times(1);
tMax = lag_times(end);

if nBasis == 1
    centers = (tMin + tMax) / 2;
    halfWidth = max([abs(tMax - tMin) / 2, lagStep]);
else
    centers = linspace(tMin, tMax, double(nBasis));
    spacing = centers(2) - centers(1);
    if spacing <= 0
        spacing = lagStep;
    end
    halfWidth = max(spacing * overlap / 2, lagStep);
end

centers = centers(:);
halfWidth = double(halfWidth);

basisMatrix = zeros(nLags, double(nBasis));
halfWidths = repmat(halfWidth, double(nBasis), 1);

for kk = 1:double(nBasis)
    span = halfWidths(kk);
    if span <= 0
        phi = zeros(nLags, 1);
    else
        phi = (lag_times - centers(kk)) / span;
    end
    mask = abs(phi) <= 1;
    if any(mask)
        tmp = zeros(nLags, 1);
        tmp(mask) = 0.5 * (cos(pi * phi(mask)) + 1);
        basisMatrix(:, kk) = tmp;
    end
end

switch lower(cfg.normalize)
    case {'l1', 'sum'}
        colSums = sum(basisMatrix, 1);
        zeroCols = colSums == 0;
        colSums(zeroCols) = 1;
        basisMatrix = basisMatrix ./ colSums;
    case {'none', ''}
        % leave as-is
    otherwise
        error('build_raised_cosine_basis:InvalidNormalize', 'normalize option ''%s'' not recognised.', cfg.normalize);
end

meta = struct();
meta.mode = 'raised_cosine';
meta.matrix = basisMatrix;
meta.n_basis = double(nBasis);
meta.centers_s = centers;
meta.half_width_s = halfWidths;
meta.normalize = cfg.normalize;
meta.lag_step_s = lagStep;
end
