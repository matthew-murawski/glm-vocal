function streams = build_streams(ev, stim)
% build binary heard/produced streams aligned to the time grid.
%
% streams = build_streams(ev, stim) returns logical column vectors marking
% bins that overlap any perceived or produced event interval on the
% stimulus timeline.

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
    streams = struct('heard_any', false(0, 1), 'produced_any', false(0, 1));
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
producedStream = false(nBins, 1);
if isempty(ev)
    streams = struct('heard_any', heardStream, 'produced_any', producedStream);
    return
end
if ~isstruct(ev) || ~all(isfield(ev, {'kind', 't_on', 't_off'}))
    error('glm:InvalidInput', 'Event structs must include kind, t_on, and t_off fields.');
end

kindCells = cellfun(@char, {ev.kind}, 'UniformOutput', false);
heardMask = strcmp(kindCells, 'perceived');
producedMask = strcmp(kindCells, 'produced');

%% rasterize each event subset
% we mark both perceived and produced call onsets with impulses so downstream kernels align to event starts.
if any(heardMask)
    heardStream = mark_event_onsets(ev(heardMask), gridStart, dt, nBins);
end
if any(producedMask)
    producedStream = mark_event_onsets(ev(producedMask), gridStart, dt, nBins);
end

%% package the result
% we expose the logical streams with a consistent struct interface used downstream.
streams = struct('heard_any', heardStream, 'produced_any', producedStream);
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
