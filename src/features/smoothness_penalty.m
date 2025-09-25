function [D, Dmap] = smoothness_penalty(colmap, cfg, lambdaCfg)
% section inputs and defaults
% we normalise configuration inputs, gather column metadata, and prepare optional lambda scaling for block penalties.
if nargin < 2 || isempty(cfg)
    cfg = struct();
end
if nargin < 3
    lambdaCfg = 1;
end
if ~isstruct(colmap)
    error('smoothness_penalty:InvalidColmap', 'colmap must be a struct describing design matrix blocks.');
end

allCols = collectAllColumns(colmap);
if isempty(allCols)
    D = sparse(0, 0);
    Dmap = struct();
    return
end
totalCols = max(allCols);

% section iterate penalised blocks
% we build second-difference stencils for each penalised kernel block while skipping scalar columns. lambda scaling weights each block's rows.
rowCounter = 0;
rowIdx = [];
colIdx = [];
valIdx = [];
Dmap = struct();

fields = fieldnames(colmap);
for ff = 1:numel(fields)
    fieldName = fields{ff};
    entry = colmap.(fieldName);
    if ~isstruct(entry) || ~isfield(entry, 'cols')
        continue
    end
    if ~shouldPenalize(fieldName, cfg)
        continue
    end

    blockCols = double(entry.cols(:)');
    nBlockCols = numel(blockCols);
    if nBlockCols < 3
        Dmap.(fieldName) = struct('rows', zeros(0, 1), 'cols', blockCols(:));
        continue
    end

    lambdaScale = resolveLambdaScale(fieldName, lambdaCfg);
    if lambdaScale == 0
        Dmap.(fieldName) = struct('rows', zeros(0, 1), 'cols', blockCols(:));
        continue
    end

    nRowsBlock = nBlockCols - 2;
    localRows = (1:nRowsBlock)';
    globalRows = rowCounter + localRows;

    colTriplets = [blockCols(1:end-2); blockCols(2:end-1); blockCols(3:end)];
    colTriplets = colTriplets(:);
    rowRep = repelem(globalRows, 3);
    stencilVals = repmat([1; -2; 1], nRowsBlock, 1) * lambdaScale;

    rowIdx = [rowIdx; rowRep]; %#ok<AGROW>
    colIdx = [colIdx; colTriplets]; %#ok<AGROW>
    valIdx = [valIdx; stencilVals]; %#ok<AGROW>

    Dmap.(fieldName) = struct('rows', globalRows, 'cols', blockCols(:), 'lambda_scale', lambdaScale);
    rowCounter = rowCounter + nRowsBlock;
end

% section sparse assembly
% after walking eligible blocks we assemble the penalty operator as a sparse matrix covering all design columns.
if rowCounter == 0
    D = sparse(0, totalCols);
else
    D = sparse(rowIdx, colIdx, valIdx, rowCounter, totalCols);
end
end

function cols = collectAllColumns(colmap)
% section column collector
% gather every column index referenced in the column map so we can size the penalty operator correctly.
cols = [];
fields = fieldnames(colmap);
for ff = 1:numel(fields)
    entry = colmap.(fields{ff});
    if ~isstruct(entry)
        continue
    end
    if numel(entry) > 1
        for ee = 1:numel(entry)
            cols = appendCols(cols, entry(ee));
        end
    else
        cols = appendCols(cols, entry);
    end
end
end

function cols = appendCols(cols, entry)
% section column helper
% append any column indices from the provided entry struct when available.
if isfield(entry, 'cols') && ~isempty(entry.cols)
    cols = [cols, double(entry.cols(:)')]; %#ok<AGROW>
end
end

function scale = resolveLambdaScale(fieldName, lambdaCfg)
% section lambda scaling
% translate scalar or struct lambda specifications into per-block square-root scaling factors.
if isnumeric(lambdaCfg)
    if ~isscalar(lambdaCfg)
        error('smoothness_penalty:LambdaScalar', 'numeric lambda must be a scalar.');
    end
    scale = scalarToScale(lambdaCfg);
    return
end

if ~isstruct(lambdaCfg)
    error('smoothness_penalty:LambdaType', 'lambda must be a scalar or struct.');
end

key = lambdaKeyForField(fieldName);
if isempty(key)
    scale = 1;
    return
end

if isfield(lambdaCfg, key)
    scale = scalarToScale(lambdaCfg.(key));
elseif isfield(lambdaCfg, fieldName)
    scale = scalarToScale(lambdaCfg.(fieldName));
else
    scale = 1;
end
end

function scale = scalarToScale(val)
% section lambda guard
% validate scalar lambda entries and convert them to square-root weights for row scaling.
val = double(val);
if ~isscalar(val) || ~isfinite(val) || val < 0
    error('smoothness_penalty:LambdaValue', 'lambda values must be finite, non-negative scalars.');
end
scale = sqrt(val);
end

function key = lambdaKeyForField(fieldName)
% section lambda mapping
% map column-map field names to canonical lambda keys.
switch fieldName
    case 'heard_any'
        key = 'heard';
    case 'produced_any'
        key = 'produced';
    case 'spike_history'
        key = 'history';
    otherwise
        key = '';
end
end

function tf = shouldPenalize(fieldName, cfg)
% section penalty selector
% decide whether a given block should receive a smoothness penalty, defaulting to kernel blocks when the config is silent.
defaultBlocks = {'heard_any', 'produced_any', 'spike_history'};
if isstruct(cfg) && isfield(cfg, fieldName)
    tf = logical(cfg.(fieldName));
else
    tf = ismember(fieldName, defaultBlocks);
end
end
