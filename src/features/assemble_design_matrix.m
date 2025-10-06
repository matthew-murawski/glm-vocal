function Xd = assemble_design_matrix(streams, states, sps, cfg, stim)
% section input validation
% confirm the stimulus timeline, streams, and configuration expose the required fields before building kernels.
mustHaveField(stim, 't');
mustHaveField(stim, 'dt');
t = stim.t(:);
dt = stim.dt;
if isempty(t)
    Xd = struct('X', sparse(0, 0), 'y', zeros(0, 1), 'colmap', struct());
    return
end
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('assemble_design_matrix:InvalidDt', 'stim.dt must be a positive finite scalar.');
end
nT = numel(t);

sps = double(sps(:));
if numel(sps) ~= nT
    error('assemble_design_matrix:SpikeLength', 'spike raster length (%d) must match stim.t (%d).', numel(sps), nT);
end

mustHaveField(streams, 'heard_any');
if numel(streams.heard_any(:)) ~= nT
    error('assemble_design_matrix:StreamLength', 'heard_any stream must match the stimulus grid.');
end

heardFields = fetch_heard_fields(streams);
heardStreams = struct();
for ii = 1:numel(heardFields)
    fieldName = heardFields{ii};
    if ~isfield(streams, fieldName)
        error('assemble_design_matrix:MissingStream', 'streams missing heard field %s.', fieldName);
    end
    heardStreams.(fieldName) = double(streams.(fieldName)(:));
    if numel(heardStreams.(fieldName)) ~= nT
        error('assemble_design_matrix:StreamLength', 'stream %s must match the stimulus grid.', fieldName);
    end
end

producedFields = fetch_produced_fields(streams);
producedStreams = struct();
for ii = 1:numel(producedFields)
    fieldName = producedFields{ii};
    if ~isfield(streams, fieldName)
        error('assemble_design_matrix:MissingStream', 'streams missing produced field %s.', fieldName);
    end
    producedStreams.(fieldName) = double(streams.(fieldName)(:));
    if numel(producedStreams.(fieldName)) ~= nT
        error('assemble_design_matrix:StreamLength', 'stream %s must match the stimulus grid.', fieldName);
    end
end

mustHaveField(states, 'convo');
mustHaveField(states, 'spon');
stateConvo = double(states.convo(:));
stateSpon = double(states.spon(:));
if numel(stateConvo) ~= nT || numel(stateSpon) ~= nT
    error('assemble_design_matrix:StateLength', 'state flags must match the stimulus grid.');
end

heardWindow = fetchWindow(cfg, 'heard_window_s');
producedWindow = fetchWindow(cfg, 'produced_window_s');
historyWindow = fetchWindow(cfg, 'history_window_s');
heardBasisCfg = fetchBasis(cfg, 'heard_basis');
producedBasisCfg = fetchBasis(cfg, 'produced_basis');

goodMask = true(nT, 1);
if isfield(stim, 'mask') && isstruct(stim.mask) && isfield(stim.mask, 'good') && ~isempty(stim.mask.good)
    maskVec = stim.mask.good(:);
    if numel(maskVec) ~= nT
        error('assemble_design_matrix:MaskLength', 'stim.mask.good must match stim.t length.');
    end
    goodMask = logical(maskVec);
end

excludeList = determine_exclusions(cfg);
skipBlock = @(name) any(strcmpi(name, excludeList));

% section block construction
blockCells = {};
colmap = struct();
colStart = 1;

if ~skipBlock('intercept')
    interceptCol = sparse((1:nT)', 1, 1, nT, 1);
    blockCells{end+1} = interceptCol; %#ok<AGROW>
    colmap.intercept = struct('cols', colStart, 'name', 'intercept');
    colStart = colStart + 1;
end

heardFieldsIncluded = {};
for ii = 1:numel(heardFields)
    fieldName = heardFields{ii};
    if skipBlock(fieldName)
        continue
    end
    z = heardStreams.(fieldName);
    [heardBlk, heardInfo] = build_basis_block(z, stim, heardWindow, heardBasisCfg, 'causal');
    blockCells{end+1} = heardBlk; %#ok<AGROW>
    nHeard = size(heardBlk, 2);
    heardCols = colStart:(colStart + nHeard - 1);
    colmap.(fieldName) = struct('cols', heardCols, 'info', heardInfo);
    heardFieldsIncluded{end+1} = fieldName; %#ok<AGROW>
    colStart = colStart + nHeard;
end
if ~isempty(heardFieldsIncluded)
    colmap.heard_fields = heardFieldsIncluded;
end

producedFieldsIncluded = {};
for ii = 1:numel(producedFields)
    fieldName = producedFields{ii};
    if skipBlock(fieldName)
        continue
    end
    z = producedStreams.(fieldName);
    [blk, info] = build_basis_block(z, stim, producedWindow, producedBasisCfg, 'symmetric');
    blockCells{end+1} = blk; %#ok<AGROW>
    nCols = size(blk, 2);
    cols = colStart:(colStart + nCols - 1);
    colmap.(fieldName) = struct('cols', cols, 'info', info);
    producedFieldsIncluded{end+1} = fieldName; %#ok<AGROW>
    colStart = colStart + nCols;
end
if ~isempty(producedFieldsIncluded)
    colmap.produced_fields = producedFieldsIncluded;
end

if ~skipBlock('states')
    stateBlk = sparse(double([stateConvo, stateSpon]));
    blockCells{end+1} = stateBlk; %#ok<AGROW>
    stateCols = colStart:(colStart + 1);
    stateMap = struct();
    stateMap.cols = stateCols;
    stateMap.names = {'convo', 'spon'};
    stateMap.convo = stateCols(1);
    stateMap.spon = stateCols(2);
    colmap.states = stateMap;
    colStart = colStart + numel(stateCols);
end

if ~skipBlock('spike_history')
    [historyBlk, historyInfo] = build_history_block(sps, stim, historyWindow);
    blockCells{end+1} = historyBlk; %#ok<AGROW>
    nHist = size(historyBlk, 2);
    historyCols = colStart:(colStart + nHist - 1);
    colmap.spike_history = struct('cols', historyCols, 'info', historyInfo);
    colStart = colStart + nHist;
end

if isempty(blockCells)
    X = sparse(nT, 0);
else
    X = cat(2, blockCells{:});
end

% section mask application
y = sps;
if any(~goodMask)
    X = X(goodMask, :);
    y = y(goodMask);
end

% section outputs
Xd = struct();
Xd.X = X;
Xd.y = y;
Xd.colmap = colmap;
end

function producedFields = fetch_produced_fields(streams)
if isstruct(streams) && isfield(streams, 'produced_fields') && ~isempty(streams.produced_fields)
    producedFields = cellstr(streams.produced_fields(:));
else
    defaults = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
    producedFields = defaults(isfield(streams, defaults));
end
end

function heardFields = fetch_heard_fields(streams)
if isstruct(streams) && isfield(streams, 'heard_fields') && ~isempty(streams.heard_fields)
    heardFields = cellstr(streams.heard_fields(:));
else
    heardFields = {'heard_any'};
end
end

function window = fetchWindow(cfg, fieldName)
if ~isstruct(cfg) || ~isfield(cfg, fieldName)
    error('assemble_design_matrix:MissingWindow', 'config missing field %s.', fieldName);
end
window = double(cfg.(fieldName));
if ~isnumeric(window) || numel(window) ~= 2
    error('assemble_design_matrix:WindowShape', 'config field %s must be a two-element numeric vector.', fieldName);
end
end

function basisCfg = fetchBasis(cfg, fieldName)
if ~isstruct(cfg) || ~isfield(cfg, fieldName)
    basisCfg = struct('kind', 'raised_cosine');
    return
end

basisCfg = cfg.(fieldName);
if ~isstruct(basisCfg)
    error('assemble_design_matrix:BasisShape', 'config field %s must be a struct.', fieldName);
end
if ~isfield(basisCfg, 'kind') || isempty(basisCfg.kind)
    basisCfg.kind = 'raised_cosine';
end
end

function mustHaveField(s, fname)
if ~isstruct(s) || ~isfield(s, fname)
    error('assemble_design_matrix:MissingField', 'struct is missing required field %s.', fname);
end
end

function excludeList = determine_exclusions(cfg)
excludeList = {};
if ~isstruct(cfg) || ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    return
end

raw = cfg.exclude_predictors;
if ischar(raw)
    raw = {raw};
elseif isstring(raw)
    raw = cellstr(raw(:));
elseif iscell(raw)
    tmp = cellfun(@(x) convertToChar(x), raw(:), 'UniformOutput', false);
    raw = tmp;
else
    error('assemble_design_matrix:ExcludeType', 'exclude_predictors must be a char vector, string, or cell array of strings.');
end

raw = raw(:);
excludeList = cellfun(@(x) lower(strtrim(x)), raw, 'UniformOutput', false);
excludeList = excludeList(~cellfun('isempty', excludeList));
if ~isempty(excludeList)
    excludeList = unique(excludeList, 'stable');
    if any(strcmp(excludeList, 'heard_any'))
        extra = {'heard_addressed', 'heard_overheard'};
        excludeList = unique([excludeList(:); extra(:)], 'stable');
    end
end
end

function out = convertToChar(val)
if isstring(val)
    out = char(val);
elseif ischar(val)
    out = val;
else
    error('assemble_design_matrix:ExcludeType', 'exclude_predictors entries must be strings.');
end
end
