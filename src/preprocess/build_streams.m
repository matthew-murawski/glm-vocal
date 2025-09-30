function streams = build_streams(ev, stim)
% build binary heard/produced streams aligned to the time grid.
%
% streams = build_streams(ev, stim) returns logical column vectors marking
% bins containing heard call onsets along with three produced-call context
% streams derived from a five second lookback window.

%% validate the stimulus grid
% we ensure the stimulus struct exposes a usable time vector and bin width before rasterizing events.
if nargin < 2
    error('glm:InvalidInput', 'Events and stimulus struct are required.');
end
if ~isstruct(stim) || ~isfield(stim, 't') || ~isfield(stim, 'dt')
    error('glm:InvalidInput', 'Stimulus struct must contain fields t and dt.');
end

% coerce the time vector into a column and capture basic grid properties
t = stim.t(:);
dt = stim.dt;
if isempty(t)
    streams = emptyStreams();
    return
end
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('glm:InvalidInput', 'Stimulus dt must be a positive finite scalar.');
end
nBins = numel(t);
gridStart = t(1);

%% partition events by kind
% we split the input events into perceived/produced lists and preallocate logical outputs.
heardStream = false(nBins, 1);
producedSpontStream = false(nBins, 1);
producedAfterHeardStream = false(nBins, 1);
producedAfterProducedStream = false(nBins, 1);
if isempty(ev)
    streams = packageStreams(heardStream, producedSpontStream, producedAfterHeardStream, producedAfterProducedStream);
    return
end
if ~isstruct(ev) || ~all(isfield(ev, {'kind', 't_on', 't_off'}))
    error('glm:InvalidInput', 'Event structs must include kind, t_on, and t_off fields.');
end

kindCells = cellfun(@char, {ev.kind}, 'UniformOutput', false);
heardMask = strcmp(kindCells, 'perceived');
producedMask = strcmp(kindCells, 'produced');

lookbackWindowS = 5.0;

%% rasterize each event subset
% we mark both perceived and produced call onsets with impulses so downstream kernels align to event starts.
if any(heardMask)
    heardStream = mark_event_onsets(ev(heardMask), gridStart, dt, nBins);
end
if any(producedMask)
    categoryIdx = classify_produced_events(ev, producedMask, lookbackWindowS);
    if ~isempty(categoryIdx.produced_spontaneous)
        producedSpontStream = mark_event_onsets(ev(categoryIdx.produced_spontaneous), gridStart, dt, nBins);
    end
    if ~isempty(categoryIdx.produced_after_heard)
        producedAfterHeardStream = mark_event_onsets(ev(categoryIdx.produced_after_heard), gridStart, dt, nBins);
    end
    if ~isempty(categoryIdx.produced_after_produced)
        producedAfterProducedStream = mark_event_onsets(ev(categoryIdx.produced_after_produced), gridStart, dt, nBins);
    end
end

%% package the result
% we expose the logical streams with a consistent struct interface used downstream.
streams = packageStreams(heardStream, producedSpontStream, producedAfterHeardStream, producedAfterProducedStream);
end

function streams = emptyStreams()
% construct an empty streams struct with all produced categories represented.
heard = false(0, 1);
produced = false(0, 1);
streams = struct( ...
    'heard_any', heard, ...
    'produced_spontaneous', produced, ...
    'produced_after_heard', produced, ...
    'produced_after_produced', produced, ...
    'produced_any', produced);
end

function streams = packageStreams(heardStream, producedSpont, producedAfterHeard, producedAfterProduced)
% bundle the individual streams while exposing an aggregate produced field for compatibility.
producedAny = producedSpont | producedAfterHeard | producedAfterProduced;
streams = struct();
streams.heard_any = heardStream;
streams.produced_spontaneous = producedSpont;
streams.produced_after_heard = producedAfterHeard;
streams.produced_after_produced = producedAfterProduced;
streams.produced_any = producedAny;
end

function stream = rasterize_events(events, gridStart, dt, nBins)
% convert an event list into a logical stream aligned to the stimulus grid.

%% initialize the stream
% we start with a false vector and exit early when no events are provided.
stream = false(nBins, 1);
if isempty(events)
    return
end

%% iterate over events
% we translate each interval into covered bin indices using inclusive start and exclusive end logic.
for ii = 1:numel(events)
    tOn = double(events(ii).t_on);
    tOff = double(events(ii).t_off);

    % skip degenerate or entirely out-of-range intervals
    if ~(isfinite(tOn) && isfinite(tOff)) || tOff <= tOn
        continue
    end

    startIdx = floor((tOn - gridStart) / dt) + 1;
    edgeTol = max([eps(dt), eps(tOff), eps(tOff - gridStart)]);
    endEdge = (tOff - gridStart) - edgeTol;
    endIdx = ceil(endEdge / dt);

    if endIdx < startIdx
        continue
    end

    startIdx = max(startIdx, 1);
    endIdx = min(endIdx, nBins);
    if startIdx > nBins || endIdx < 1
        continue
    end

    stream(startIdx:endIdx) = true;
end
end

function impulses = mark_event_onsets(events, gridStart, dt, nBins)
% convert produced call onsets into a sparse impulse train aligned to the grid.

%% initialize the impulse stream
% we start with zeros so only bins containing call onsets receive a one.
impulses = false(nBins, 1);
if isempty(events)
    return
end

%% iterate over events
% we locate the bin containing each onset and flip it on, clipping out-of-range onsets.
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
