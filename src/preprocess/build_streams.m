function streams = build_streams(ev, stim, cfg)
% build binary heard/produced streams aligned to the time grid.
%
% streams = build_streams(ev, stim, cfg) returns logical column vectors marking
% heard call onsets and produced-call categories. the produced categories are
% determined by cfg.produced_split_mode, which can be either 'context' (the
% legacy behaviour) or 'call_type'.

%% validate the stimulus grid
if nargin < 2
    error('glm:InvalidInput', 'Events and stimulus struct are required.');
end
if nargin < 3
    cfg = struct();
end
if ~isstruct(stim) || ~isfield(stim, 't') || ~isfield(stim, 'dt')
    error('glm:InvalidInput', 'Stimulus struct must contain fields t and dt.');
end

t = stim.t(:);
dt = stim.dt;
if isempty(t)
    streams = empty_streams(cfg, 0);
    return
end
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('glm:InvalidInput', 'Stimulus dt must be a positive finite scalar.');
end

nBins = numel(t);
gridStart = t(1);

%% partition events by kind
heardStream = false(nBins, 1);
if isempty(ev)
    streams = empty_streams(cfg, nBins);
    return
end
if ~isstruct(ev) || ~all(isfield(ev, {'kind', 't_on', 't_off'}))
    error('glm:InvalidInput', 'Event structs must include kind, t_on, and t_off fields.');
end

ev = ev(:);
kindCells = cellfun(@char, {ev.kind}, 'UniformOutput', false);
heardMask = strcmp(kindCells, 'perceived');
producedMask = strcmp(kindCells, 'produced');

lookbackWindowS = 5.0;

%% rasterize heard events
if any(heardMask)
    heardStream = mark_event_onsets(ev(heardMask), gridStart, dt, nBins);
end

%% classify produced events according to the requested mode
categoryIdx = classify_produced_events(ev, producedMask, lookbackWindowS, cfg);
producedNames = categoryIdx.names;
if isempty(producedNames)
    producedNames = {}; %#ok<NASGU>
end

producedStreams = struct();
producedAny = false(nBins, 1);
for ii = 1:numel(producedNames)
    fieldName = producedNames{ii};
    idx = categoryIdx.(fieldName);
    if isempty(idx)
        producedStream = false(nBins, 1);
    else
        producedStream = mark_event_onsets(ev(idx), gridStart, dt, nBins);
    end
    producedStreams.(fieldName) = producedStream;
    producedAny = producedAny | producedStream;
end

%% package outputs
streams = package_streams(heardStream, producedStreams, producedNames, producedAny, nBins);
end

function streams = empty_streams(cfg, nBins)
modeCategories = classify_produced_events(repmat(empty_event_record(), 0, 1), false(0, 1), 5.0, cfg);
producedNames = modeCategories.names;
heard = false(nBins, 1);
producedAny = false(nBins, 1);
producedStreams = struct();
for ii = 1:numel(producedNames)
    producedStreams.(producedNames{ii}) = false(nBins, 1);
end
streams = package_streams(heard, producedStreams, producedNames, producedAny, nBins);
end

function streams = package_streams(heardStream, producedStreams, producedNames, producedAny, nBins)
streams = struct();
streams.heard_any = heardStream;
for ii = 1:numel(producedNames)
    fieldName = producedNames{ii};
    streams.(fieldName) = producedStreams.(fieldName);
end
streams.produced_any = producedAny;
streams.produced_fields = producedNames(:)';

if ~isfield(streams, 'produced_any')
    streams.produced_any = false(nBins, 1);
end
end

function impulses = mark_event_onsets(events, gridStart, dt, nBins)
impulses = false(nBins, 1);
if isempty(events)
    return
end

events = events(:);
for ii = 1:numel(events)
    tOn = double(events(ii).t_on);
    if ~isfinite(tOn)
        continue
    end

    if isfield(events, 't_off')
        tOff = double(events(ii).t_off);
        if isfinite(tOff) && tOff <= tOn
            continue
        end
    end

    offset = (tOn - gridStart) / dt;
    tol = max([1e-9, eps(offset), eps(dt)]);
    idx = floor(offset + tol) + 1;
    if idx < 1 || idx > nBins
        continue
    end

    impulses(idx) = true;
end
end

function rec = empty_event_record()
rec = struct('kind', '', 't_on', [], 't_off', [], 'label', "");
end
