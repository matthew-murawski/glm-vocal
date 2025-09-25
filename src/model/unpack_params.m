function kernels = unpack_params(wmap, colmap, cfg, stim)
% section inputs and setup
% unpack the map weights into named kernel blocks using the column map metadata.
unused = {cfg, stim}; %#ok<NASGU>
w = wmap.w(:);
kernels = struct();

% section intercept
% extract the scalar intercept term when present.
if isfield(colmap, 'intercept') && isfield(colmap.intercept, 'cols')
    kernels.intercept = w(colmap.intercept.cols);
end

% section heard kernel
% slice the heard stimulus kernel weights and attach lag metadata.
if isfield(colmap, 'heard_any') && isfield(colmap.heard_any, 'cols')
    kernels.heard_any = buildKernelBlock(w, colmap.heard_any);
end

% section produced kernel
% slice the produced kernel weights (symmetric window) with metadata.
if isfield(colmap, 'produced_any') && isfield(colmap.produced_any, 'cols')
    kernels.produced_any = buildKernelBlock(w, colmap.produced_any);
end

% section spike-history kernel
% slice the spike-history kernel weights and metadata if present.
if isfield(colmap, 'spike_history') && isfield(colmap.spike_history, 'cols')
    kernels.spike_history = buildKernelBlock(w, colmap.spike_history);
end

% section state weights
% expose state weights with human-readable names when available.
if isfield(colmap, 'states') && isfield(colmap.states, 'cols')
    stateCols = colmap.states.cols(:)';
    stateWeights = w(stateCols);
    stateStruct = struct();
    stateStruct.weights = stateWeights;
    if isfield(colmap.states, 'names')
        stateStruct.names = colmap.states.names;
        for ii = 1:numel(colmap.states.names)
            name = colmap.states.names{ii};
            stateStruct.(name) = stateWeights(ii);
        end
    end
    kernels.states = stateStruct;
end
end

function block = buildKernelBlock(w, blockDef)
% section kernel helper
% create a kernel struct populated with weights and any metadata supplied in colmap.
weights = w(blockDef.cols(:));
block = struct();
block.weights = weights;
if isfield(blockDef, 'info') && isstruct(blockDef.info)
    info = blockDef.info;
    infoFields = fieldnames(info);
    for ii = 1:numel(infoFields)
        fieldName = infoFields{ii};
        block.(fieldName) = info.(fieldName);
    end
end
end
