function orchestrate_eta_sanity_check(animalID, sessionNumber, dataPath, ChanGeo, etaWindow_s, target_dt)
% orchestrate_eta_sanity_check coordinate loading, eta computation, and plotting.
% this function loops over event types and channel grid positions, loads
% the appropriate high-gamma trace, computes an event-triggered average
% around onsets for the selected event type, and plots results in a grid.

%% section: validate inputs and derive grid
% we keep validation light here and rely on downstream loaders to error if inputs are bad.
if nargin < 6
    error('orchestrate_eta_sanity_check:InvalidInput', 'expected 6 inputs: animalID, sessionNumber, dataPath, ChanGeo, etaWindow_s, target_dt');
end

if ~ismatrix(ChanGeo)
    error('orchestrate_eta_sanity_check:InvalidInput', 'ChanGeo must be a 2d matrix of channel ids');
end
[nRows, nCols] = size(ChanGeo);

% resolve hg path once; optionally preload simple layout for speed
hgPath = resolve_hg_path(dataPath);
[simpleLFPMat, simpleTimeVec] = try_preload_simple_hg(hgPath, target_dt);

% optionally start a parallel pool if available
ensure_parallel_pool();

%% section: locate and load labels once (heard + produced)
% prefer files under P.repo_root/data/Label Files/S<session> with
% patterns <animal>_S<session>_heard_*.txt and <animal>_S<session>_produced_*.txt.
events = load_session_events(dataPath, animalID, sessionNumber);
eventKinds = cellfun(@char, {events.kind}, 'UniformOutput', false);

% build onsets for produced and perceived, then split perceived into addressed/overheard
producedOnsets = [];
perceivedOnsets = [];
if ~isempty(events)
    if any(strcmpi(eventKinds, 'produced'))
        producedOnsets = [events(strcmpi(eventKinds, 'produced')).t_on];
    end
    if any(strcmpi(eventKinds, 'perceived'))
        perceivedOnsets = [events(strcmpi(eventKinds, 'perceived')).t_on];
    end
end
producedOnsets = producedOnsets(:);
perceivedOnsets = perceivedOnsets(:);

addressWindow_s = 5.0;
addressedOnsets = [];
overheardOnsets = [];
if isempty(perceivedOnsets)
    addressedOnsets = [];
    overheardOnsets = [];
elseif isempty(producedOnsets)
    addressedOnsets = [];
    overheardOnsets = perceivedOnsets;
else
    % classify each perceived onset by proximity to a future produced onset
    Np = numel(perceivedOnsets);
    isAddressed = false(Np,1);
    for ii = 1:Np
        t = perceivedOnsets(ii);
        dt = producedOnsets - t;
        dtFuture = dt(dt >= 0);
        if ~isempty(dtFuture) && min(dtFuture) <= addressWindow_s
            isAddressed(ii) = true;
        end
    end
    addressedOnsets = perceivedOnsets(isAddressed);
    overheardOnsets = perceivedOnsets(~isAddressed);
end

% define event types to process and provide onset lookup
eventTypes = {'produced', 'addressed', 'overheard'};
onsetsByType = struct('produced', producedOnsets, 'addressed', addressedOnsets, 'overheard', overheardOnsets);

%% section: iterate event types and build figures
% for each event type, create a figure with a grid of subplots matching the channel geometry.
for eIdx = 1:numel(eventTypes)
    eventType = eventTypes{eIdx};

    % collect onsets for the current event type
    if isfield(onsetsByType, eventType)
        event_onsets_list = onsetsByType.(eventType);
    else
        event_onsets_list = [];
    end

    % create figure and set a descriptive name
    figTitle = sprintf('ETA - %s Calls', capitalize_first(eventType));
    f = figure('Name', figTitle); %#ok<NASGU>

    % prebuild index maps for valid channels
    validMask = isfinite(ChanGeo) & ChanGeo > 0;
    validSubIdx = find(validMask(:));
    validChans = ChanGeo(validMask);
    nValid = numel(validChans);
    eta_time_cell = cell(nValid, 1);
    eta_trace_cell = cell(nValid, 1);

    % compute eta per-channel in parallel when possible
    parfor k = 1:nValid
        chanID = validChans(k);

        % load channel power/time either from preloaded matrix or via loader
        if ~isempty(simpleLFPMat)
            colIdx = channel_index(simpleLFPMat, chanID);
            power = double(simpleLFPMat(:, colIdx));
            tvec = simpleTimeVec;
        else
            hgLocal = [];
            try
                hgLocal = load_high_gamma(hgPath, 'channel_id', chanID);
            catch ME
                simpleHg = try_load_hg_simple(hgPath, chanID, target_dt);
                if ~isempty(simpleHg)
                    hgLocal = simpleHg;
                else
                    rethrow(ME);
                end
            end
            power = hgLocal.power;
            tvec = hgLocal.t;
        end

        [eta_tr, eta_tm] = calculate_eta_from_continuous(power, tvec, event_onsets_list, etaWindow_s, target_dt);
        eta_time_cell{k} = eta_tm;
        eta_trace_cell{k} = eta_tr;
    end

    % map from linear subplot index to computed cell index
    gridIndexMap = zeros(nRows * nCols, 1);
    for k = 1:nValid
        gridIndexMap(validSubIdx(k)) = k;
    end

    % create all subplots; blank out invalid, plot computed channels
    for r = 1:nRows
        for c = 1:nCols
            subIdx = (r - 1) * nCols + c;
            chanID = ChanGeo(r, c);
            ax = subplot(nRows, nCols, subIdx); %#ok<LAXES>

            if ~isfinite(chanID) || chanID <= 0
                axis(ax, 'off');
                title(ax, '');
                continue
            end

            k = gridIndexMap(subIdx);
            if k > 0
                eta_time = eta_time_cell{k};
                eta_trace = eta_trace_cell{k};
                plot_eta_trace(ax, eta_time, eta_trace);
                title(ax, sprintf('Channel %d', chanID));
            else
                axis(ax, 'off');
                title(ax, '');
            end
        end
    end

    % add a main title spanning the figure
    sgtitle(sprintf('High-Gamma ETA (%s) | %s | Session %d', capitalize_first(eventType), animalID, sessionNumber));
end

end

%% helpers (internal)
function out = capitalize_first(str)
% shallow helper to capitalize first character for nice titles
if isempty(str)
    out = str;
else
    out = [upper(str(1)), str(2:end)];
end
end

function events = load_session_events(dataPath, animalID, sessionNumber)
% assemble a combined events struct array for the session by loading heard and
% produced label files when available. if not found, attempt a best-effort load
% from any label file in the provided path.

% try to get the repo root from base workspace P without modifying it
repoRoot = try_get_repo_root();

% choose a base directory to search for labels
baseDir = '';
if ~isempty(repoRoot)
    baseDir = fullfile(repoRoot, 'data', 'Label Files', sprintf('S%d', sessionNumber));
elseif exist(dataPath, 'dir') == 7
    baseDir = dataPath;
elseif exist(dataPath, 'file') == 2
    baseDir = fileparts(dataPath);
end

% initialize containers
heardEvents = struct('kind', {}, 't_on', {}, 't_off', {}, 'label', {});
producedEvents = heardEvents;

% search for canonical filenames first
heardFile = '';
producedFile = '';
if ~isempty(baseDir) && exist(baseDir, 'dir') == 7
    prefix = sprintf('%s_S%d_', char(string(animalID)), sessionNumber);
    cand = dir(fullfile(baseDir, [prefix '*heard*.txt']));
    if ~isempty(cand), heardFile = fullfile(cand(1).folder, cand(1).name); end
    cand = dir(fullfile(baseDir, [prefix '*produced*.txt']));
    if ~isempty(cand), producedFile = fullfile(cand(1).folder, cand(1).name); end

    % relax patterns if not found
    if isempty(heardFile)
        cand = dir(fullfile(baseDir, '*heard*.txt'));
        if ~isempty(cand), heardFile = fullfile(cand(1).folder, cand(1).name); end
    end
    if isempty(producedFile)
        cand = dir(fullfile(baseDir, '*produced*.txt'));
        if ~isempty(cand), producedFile = fullfile(cand(1).folder, cand(1).name); end
    end
end

% attempt to load with explicit default kinds if files exist
try
    if ~isempty(heardFile) && exist(heardFile, 'file') == 2
        heardEvents = load_labels(heardFile, 'perceived');
    end
catch %#ok<CTCH>
    heardEvents = struct('kind', {}, 't_on', {}, 't_off', {}, 'label', {});
end

try
    if ~isempty(producedFile) && exist(producedFile, 'file') == 2
        producedEvents = load_labels(producedFile, 'produced');
    end
catch %#ok<CTCH>
    producedEvents = struct('kind', {}, 't_on', {}, 't_off', {}, 'label', {});
end

% fallback: if neither set loaded, try any label file nearby
if isempty(heardEvents) && isempty(producedEvents)
    probeDir = baseDir;
    if isempty(probeDir) && exist(dataPath, 'dir') == 7
        probeDir = dataPath;
    end
    if ~isempty(probeDir) && exist(probeDir, 'dir') == 7
        dd = dir(fullfile(probeDir, '*.txt'));
        if ~isempty(dd)
            try
                tmp = load_labels(fullfile(dd(1).folder, dd(1).name));
                heardEvents = tmp;  %#ok<NASGU>
                % let tmp carry kinds if present; otherwise, this may error which is fine
            catch %#ok<CTCH>
                % leave empty on failure
            end
        end
    end
end

% concatenate and return (may be empty)
events = [heardEvents; producedEvents];
if isempty(events)
    events = struct('kind', {}, 't_on', {}, 't_off', {}, 'label', {});
end
end

function repoRoot = try_get_repo_root()
% safely fetch P.repo_root from the base workspace without modifying P.
repoRoot = '';
try
    P = evalin('base', 'P'); %#ok<NASGU>
    repoRoot = evalin('base', 'P.repo_root');
catch %#ok<CTCH>
    repoRoot = '';
end
if isstring(repoRoot)
    repoRoot = char(repoRoot);
end
end

function hgPath = resolve_hg_path(dataPath)
% try to resolve a high-gamma mat file from dataPath. if a file is supplied, use it;
% if a directory is supplied, search for a mat file with hg-like name, else any mat.
if isstring(dataPath)
    dataPath = char(dataPath);
end

if exist(dataPath, 'file') == 2
    hgPath = dataPath;
    return
end

if exist(dataPath, 'dir') == 7
    % try common patterns for high gamma files
    patterns = {fullfile(dataPath, '*high*gamma*.mat'), fullfile(dataPath, '*hg*.mat'), fullfile(dataPath, '*lfp*.mat')};
    for k = 1:numel(patterns)
        dd = dir(patterns{k});
        if ~isempty(dd)
            hgPath = fullfile(dd(1).folder, dd(1).name);
            return
        end
    end
    % fallback: any mat file
    dd = dir(fullfile(dataPath, '*.mat'));
    if ~isempty(dd)
        hgPath = fullfile(dd(1).folder, dd(1).name);
        return
    end
end

% last resort: return input; downstream will raise a clear file-not-found error
hgPath = dataPath;
end

function hg = try_load_hg_simple(hgPath, chanID, target_dt)
% attempt to load a simpler high-gamma file layout where the MAT contains a
% numeric variable `lfp_data` of size [nSamples x nChannels] only. if present,
% construct a minimal hg struct with a synthetic time vector based on target_dt.

hg = [];
if exist(hgPath, 'file') ~= 2
    return
end

info = whos('-file', hgPath, 'lfp_data');
if isempty(info) || ~isnumeric_placeholder(info)
    return
end

S = load(hgPath, 'lfp_data');
lfp = S.lfp_data;
if ~isnumeric(lfp) || ndims(lfp) ~= 2
    return
end

[nSamples, nChannels] = size(lfp);
if chanID < 1 || chanID > nChannels
    return
end

power = double(lfp(:, chanID));
t = (0:nSamples-1)' * target_dt;  % assume sample step equals target dt
fs = 1/target_dt;

hg = struct('power', power(:), 'fs', fs, 't', t(:), 'session_id', '', 'channel_id', chanID);
end

function tf = isnumeric_placeholder(info)
% helper: determine if loaded variable is numeric by inspecting whos output
try
    tf = ismember(info.class, {'double','single','uint16','uint32','uint8','int16','int32','int8'});
catch
    tf = false;
end
end

function ensure_parallel_pool()
% start a parallel pool if a license is available and none is active.
try
    if license('test', 'Distrib_Computing_Toolbox')
        p = gcp('nocreate');
        if isempty(p)
            parpool('local'); %#ok<*NOPRT>
        end
    end
catch %#ok<CTCH>
    % ignore failures; fall back to serial computation
end
end

function [lfpMat, tvec] = try_preload_simple_hg(hgPath, target_dt)
% fast path: if MAT contains numeric lfp_data [nSamples x nChannels], preload it once.
lfpMat = [];
tvec = [];
if exist(hgPath, 'file') ~= 2
    return
end
info = whos('-file', hgPath, 'lfp_data');
if isempty(info) || strcmp(info.class, 'struct')
    return
end
S = load(hgPath, 'lfp_data');
lfp = S.lfp_data;
if ~isnumeric(lfp) || ndims(lfp) ~= 2
    return
end
nSamples = size(lfp, 1);
tvec = (0:nSamples-1)' * target_dt;
lfpMat = lfp;
end

function idx = channel_index(lfpMat, chanID)
% assume channels are 1-indexed columns; validate range
nChannels = size(lfpMat, 2);
if chanID < 1 || chanID > nChannels
    error('orchestrate_eta_sanity_check:InvalidChannel', 'channel id %d out of range [1,%d]', chanID, nChannels);
end
idx = chanID;
end
