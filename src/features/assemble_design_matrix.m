function Xd = assemble_design_matrix(streams, states, response_data, cfg, stim)
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

% determine response type from config (default to 'spikes' for backward compatibility)
if isstruct(cfg) && isfield(cfg, 'response_type')
    response_type = cfg.response_type;
else
    response_type = 'spikes';
end

% validate response_data based on response type
response_data = double(response_data(:));
if numel(response_data) ~= nT
    error('assemble_design_matrix:ResponseLength', 'response data length (%d) must match stim.t (%d).', numel(response_data), nT);
end

% for backward compatibility, still accept 'sps' as variable name
sps = response_data;

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

% fetch history window based on response type
if strcmpi(response_type, 'lfp') && isfield(cfg, 'lfp') && isfield(cfg.lfp, 'history_window_s')
    historyWindow = fetchWindow(cfg.lfp, 'history_window_s');
else
    historyWindow = fetchWindow(cfg, 'history_window_s');
end

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

% build history block based on response type
if strcmpi(response_type, 'spikes')
    history_block_name = 'spike_history';
    if ~skipBlock(history_block_name)
        [historyBlk, historyInfo] = build_history_block(response_data, stim, historyWindow);
        blockCells{end+1} = historyBlk; %#ok<AGROW>
        nHist = size(historyBlk, 2);
        historyCols = colStart:(colStart + nHist - 1);
        colmap.spike_history = struct('cols', historyCols, 'info', historyInfo);
        colStart = colStart + nHist;
    end
elseif strcmpi(response_type, 'lfp')
    history_block_name = 'lfp_history';
    if ~skipBlock(history_block_name)
        % check if LFP history should use basis projection
        use_lfp_basis = false;
        if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'history_basis') && ...
           isstruct(cfg.lfp.history_basis) && isfield(cfg.lfp.history_basis, 'kind') && ...
           ~strcmpi(cfg.lfp.history_basis.kind, 'raw')
            use_lfp_basis = true;
        end

        [historyBlk, historyInfo] = build_lfp_history_block(response_data, stim, historyWindow, use_lfp_basis);
        blockCells{end+1} = historyBlk; %#ok<AGROW>
        nHist = size(historyBlk, 2);
        historyCols = colStart:(colStart + nHist - 1);
        colmap.lfp_history = struct('cols', historyCols, 'info', historyInfo);
        colStart = colStart + nHist;
    end
else
    error('assemble_design_matrix:InvalidResponseType', 'response_type must be "spikes" or "lfp", got: %s', response_type);
end

if isempty(blockCells)
    X = sparse(nT, 0);
else
    X = cat(2, blockCells{:});
end

% section mask application
y = response_data;
if any(~goodMask)
    X = X(goodMask, :);
    y = y(goodMask);
end

% section outputs
Xd = struct();
Xd.X = X;
Xd.y = y;
Xd.colmap = colmap;
Xd.response_type = response_type;
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
