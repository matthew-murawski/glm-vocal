function Xd = assemble_design_matrix(streams, states, sps, cfg, stim)
% section input validation
% we confirm timelines, streams, and configuration structs expose the fields needed to build kernel blocks.
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
mustHaveField(streams, 'produced_spontaneous');
mustHaveField(streams, 'produced_after_heard');
mustHaveField(streams, 'produced_after_produced');

heard = double(streams.heard_any(:));
producedSpont = double(streams.produced_spontaneous(:));
producedAfterHeard = double(streams.produced_after_heard(:));
producedAfterProduced = double(streams.produced_after_produced(:));
if numel(heard) ~= nT || numel(producedSpont) ~= nT || numel(producedAfterHeard) ~= nT || numel(producedAfterProduced) ~= nT
    error('assemble_design_matrix:StreamLength', 'stream lengths must match the stimulus grid.');
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

% section block construction
% we build each regressor block using the feature helpers so we can concatenate them into the design matrix.
interceptCol = sparse((1:nT)', 1, 1, nT, 1);
[heardBlk, heardInfo] = build_basis_block(heard, stim, heardWindow, heardBasisCfg, 'causal');
[producedSpontBlk, producedSpontInfo] = build_basis_block(producedSpont, stim, producedWindow, producedBasisCfg, 'symmetric');
[producedAfterHeardBlk, producedAfterHeardInfo] = build_basis_block(producedAfterHeard, stim, producedWindow, producedBasisCfg, 'symmetric');
[producedAfterProducedBlk, producedAfterProducedInfo] = build_basis_block(producedAfterProduced, stim, producedWindow, producedBasisCfg, 'symmetric');
stateBlk = sparse(double([stateConvo, stateSpon]));
[historyBlk, historyInfo] = build_history_block(sps, stim, historyWindow);

X = [interceptCol, heardBlk, producedSpontBlk, producedAfterHeardBlk, producedAfterProducedBlk, stateBlk, historyBlk];

% section mask application
% we drop rows marked bad by the stimulus mask so downstream fitting only sees valid bins.
y = sps;
if any(~goodMask)
    X = X(goodMask, :);
    y = y(goodMask);
end

% section column mapping
% we record the column indices associated with each block along with the helper metadata for downstream unpacking.
colStart = 1;
colmap = struct();

colmap.intercept = struct('cols', colStart, 'name', 'intercept');
colStart = colStart + 1;

nHeard = size(heardBlk, 2);
heardCols = colStart:(colStart + nHeard - 1);
colmap.heard_any = struct('cols', heardCols, 'info', heardInfo);
colStart = colStart + nHeard;

nProducedSpont = size(producedSpontBlk, 2);
producedSpontCols = colStart:(colStart + nProducedSpont - 1);
colmap.produced_spontaneous = struct('cols', producedSpontCols, 'info', producedSpontInfo);
colStart = colStart + nProducedSpont;

nProducedAfterHeard = size(producedAfterHeardBlk, 2);
producedAfterHeardCols = colStart:(colStart + nProducedAfterHeard - 1);
colmap.produced_after_heard = struct('cols', producedAfterHeardCols, 'info', producedAfterHeardInfo);
colStart = colStart + nProducedAfterHeard;

nProducedAfterProduced = size(producedAfterProducedBlk, 2);
producedAfterProducedCols = colStart:(colStart + nProducedAfterProduced - 1);
colmap.produced_after_produced = struct('cols', producedAfterProducedCols, 'info', producedAfterProducedInfo);
colStart = colStart + nProducedAfterProduced;

stateCols = colStart:(colStart + 1);
stateMap = struct();
stateMap.cols = stateCols;
stateMap.names = {'convo', 'spon'};
stateMap.convo = stateCols(1);
stateMap.spon = stateCols(2);
colmap.states = stateMap;
colStart = colStart + numel(stateCols);

nHist = size(historyBlk, 2);
historyCols = colStart:(colStart + nHist - 1);
colmap.spike_history = struct('cols', historyCols, 'info', historyInfo);

% section outputs
% we package the sparse design matrix, binned spikes, and column mapping for downstream fitting and analysis.
Xd = struct();
Xd.X = X;
Xd.y = y;
Xd.colmap = colmap;
end

function window = fetchWindow(cfg, fieldName)
% section window helper
% this helper extracts a two-element window vector from the config struct or raises an informative error.
if ~isstruct(cfg) || ~isfield(cfg, fieldName)
    error('assemble_design_matrix:MissingWindow', 'config missing field %s.', fieldName);
end
window = double(cfg.(fieldName));
if ~isnumeric(window) || numel(window) ~= 2
    error('assemble_design_matrix:WindowShape', 'config field %s must be a two-element numeric vector.', fieldName);
end
end

function basisCfg = fetchBasis(cfg, fieldName)
% section basis helper
% pull a basis configuration struct or fall back to the default raised-cosine family.
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
% section struct guard
% small helper to ensure required struct fields exist before continuing computations.
if ~isstruct(s) || ~isfield(s, fname)
    error('assemble_design_matrix:MissingField', 'struct is missing required field %s.', fname);
end
end
