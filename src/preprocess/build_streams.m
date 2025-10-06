function streams = build_streams(ev, stim, cfg)
% build binary heard/produced streams aligned to the time grid.
%
% streams = build_streams(ev, stim, cfg) returns logical column vectors marking
% heard call onsets and produced-call categories. the produced categories are
% determined by cfg.produced_split_mode, which can be either 'context' (the
% legacy behaviour) or 'call_type'. addressed and overheard heard-call streams
% are provided alongside the aggregate heard_any stream.

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
addressedWindowS = resolve_window(cfg, 'perceived_addressed_window_s', 4.0);
silenceWindowS = resolve_window(cfg, 'perceived_overheard_silence_s', 5.0);

%% rasterize heard events with addressed/overheard splits
heardCategories = classify_heard_events(ev, heardMask, producedMask, addressedWindowS, silenceWindowS);
heardStreams = struct();
heardStreams.heard_any = false(nBins, 1);
heardStreams.heard_addressed = false(nBins, 1);
heardStreams.heard_overheard = false(nBins, 1);

if any(heardMask)
    heardStreams.heard_any = mark_event_onsets(ev(heardMask), gridStart, dt, nBins);
end
if ~isempty(heardCategories.addressed)
    heardStreams.heard_addressed = mark_event_onsets(ev(heardCategories.addressed), gridStart, dt, nBins);
end
if ~isempty(heardCategories.overheard)
    heardStreams.heard_overheard = mark_event_onsets(ev(heardCategories.overheard), gridStart, dt, nBins);
end

heardStreams.heard_fields = {'heard_addressed', 'heard_overheard'};

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
streams = package_streams(heardStreams, producedStreams, producedNames, producedAny, nBins);
end

function streams = empty_streams(cfg, nBins)
modeCategories = classify_produced_events(repmat(empty_event_record(), 0, 1), false(0, 1), 5.0, cfg);
producedNames = modeCategories.names;
heardStreams = struct();
heardStreams.heard_any = false(nBins, 1);
heardStreams.heard_addressed = false(nBins, 1);
heardStreams.heard_overheard = false(nBins, 1);
heardStreams.heard_fields = {'heard_addressed', 'heard_overheard'};
producedAny = false(nBins, 1);
producedStreams = struct();
for ii = 1:numel(producedNames)
    producedStreams.(producedNames{ii}) = false(nBins, 1);
end
streams = package_streams(heardStreams, producedStreams, producedNames, producedAny, nBins);
end

function streams = package_streams(heardStreams, producedStreams, producedNames, producedAny, nBins)
streams = struct();
streams.heard_any = heardStreams.heard_any;
streams.heard_addressed = heardStreams.heard_addressed;
streams.heard_overheard = heardStreams.heard_overheard;
streams.heard_fields = heardStreams.heard_fields;
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

function window = resolve_window(cfg, fieldName, defaultVal)
% section resolve window
% pick the configured scalar duration when available, otherwise fall back to the provided default.
window = defaultVal;
if ~isstruct(cfg) || ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    return
end

candidate = double(cfg.(fieldName));
if ~isscalar(candidate) || ~isfinite(candidate) || candidate < 0
    error('build_streams:InvalidWindow', 'config field %s must be a non-negative scalar.', fieldName);
end
window = candidate;
end

function heardCategories = classify_heard_events(events, heardMask, producedMask, addressedWindowS, silenceWindowS)
% section heard classification
% assign perceived events to addressed or overheard bins using produced-call proximity rules.
heardCategories = struct('addressed', zeros(0, 1), 'overheard', zeros(0, 1));
if isempty(events) || ~any(heardMask)
    return
end

heardIdx = find(logical(heardMask(:)));
producedIdx = find(logical(producedMask(:)));

if isempty(producedIdx)
    heardCategories.overheard = heardIdx;
    return
end

starts = arrayfun(@(idx) double(events(idx).t_on), producedIdx);
stops = arrayfun(@(idx) lookup_stop(events(idx)), producedIdx);

validMask = isfinite(starts);
starts = starts(validMask);
stops = stops(validMask);
if isempty(starts)
    heardCategories.overheard = heardIdx;
    return
end

[starts, stops] = adjust_invalid_stops(starts, stops);

addressedList = [];
overheardList = [];

for ii = 1:numel(heardIdx)
    evIdx = heardIdx(ii);
    tOn = double(events(evIdx).t_on);
    if ~isfinite(tOn)
        continue
    end

    distance = min_distance_to_intervals(tOn, starts, stops);
    if distance <= addressedWindowS
        addressedList(end+1, 1) = evIdx; %#ok<AGROW>
    elseif distance >= silenceWindowS
        overheardList(end+1, 1) = evIdx; %#ok<AGROW>
    end
end

heardCategories.addressed = addressedList;
heardCategories.overheard = overheardList;
end

function stopVal = lookup_stop(event)
stopVal = double(event.t_off);
if ~isfinite(stopVal)
    stopVal = double(event.t_on);
end
end

function [startsOut, stopsOut] = adjust_invalid_stops(startsIn, stopsIn)
startsOut = startsIn(:);
stopsOut = stopsIn(:);
for ii = 1:numel(stopsOut)
    if ~isfinite(stopsOut(ii))
        stopsOut(ii) = startsOut(ii);
    end
    if stopsOut(ii) < startsOut(ii)
        tmp = stopsOut(ii);
        stopsOut(ii) = startsOut(ii);
        startsOut(ii) = tmp;
    end
end
end

function distance = min_distance_to_intervals(tPoint, starts, stops)
% section distance helper
% compute the minimal absolute distance between a time point and any produced-call interval.
if isempty(starts)
    distance = inf;
    return
end

distance = inf;
for jj = 1:numel(starts)
    s = starts(jj);
    e = stops(jj);
    if ~isfinite(s)
        continue
    end
    if ~isfinite(e)
        e = s;
    end
    if e < s
        tmp = e;
        e = s;
        s = tmp;
    end

    if tPoint >= s && tPoint <= e
        distance = 0;
        return
    elseif tPoint < s
        distance = min(distance, s - tPoint);
    else
        distance = min(distance, tPoint - e);
    end
end
end
