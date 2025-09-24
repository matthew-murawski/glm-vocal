function validate_inputs(sp, ev, sessionDuration, buffer)
%VALIDATE_INPUTS Validate spikes and events prior to processing.
%   VALIDATE_INPUTS(sp, ev, sessionDuration, buffer) ensures that spike
%   times are finite and non-decreasing, events are well-formed, and all
%   times fall within an acceptable padded session window.

% normalize optional arguments for later comparisons
if nargin < 4 || isempty(buffer)
    buffer = 0;
end
if nargin < 3 || isempty(sessionDuration)
    sessionDuration = inf;
end

% check that spike data is present and conforms to expectations
validate_spikes(sp);

% iterate over events and confirm their timing data
validate_events(ev);

% enforce global time bounds if a duration was provided
validate_time_bounds(sp, ev, sessionDuration, buffer);
end

function validate_spikes(sp)
% ensure spike struct exists and contains finite, sorted times
if ~isstruct(sp) || ~isfield(sp, 'spike_times')
    error('glm:InvalidSpikeTimes', 'Spike data must include a spike_times vector.');
end
spikeTimes = sp.spike_times;
if isempty(spikeTimes)
    return
end
if ~isa(spikeTimes, 'double') || ~isvector(spikeTimes) || any(~isfinite(spikeTimes))
    error('glm:InvalidSpikeTimes', 'Spike times must be finite doubles.');
end
spikeTimes = spikeTimes(:);
if any(diff(spikeTimes) < 0)
    error('glm:InvalidSpikeTimes', 'Spike times must be non-decreasing.');
end
end

function validate_events(ev)
% accept empty events and bypass further checks
if isempty(ev)
    return
end

% require struct array (not table/cell) with time fields present
if ~isstruct(ev)
    error('glm:InvalidLabels', 'Events must be provided as a struct array.');
end
requiredFields = {'t_on', 't_off'};
for idx = 1:numel(requiredFields)
    if ~all(isfield(ev, requiredFields{idx}))
        error('glm:InvalidLabels', 'Event struct array must include %s.', requiredFields{idx});
    end
end

% verify times are finite and respect ordering for each event
for ii = 1:numel(ev)
    tOn = ev(ii).t_on;
    tOff = ev(ii).t_off;
    if ~isscalar(tOn) || ~isscalar(tOff) || ~isfinite(tOn) || ~isfinite(tOff)
        error('glm:InvalidLabels', 'Event times must be finite scalars.');
    end
    if tOff < tOn
        error('glm:InvalidLabels', 'Event end time must be greater than or equal to start time.');
    end
end
end

function validate_time_bounds(sp, ev, sessionDuration, buffer)
% skip checks when no finite bound exists
if ~isfinite(sessionDuration)
    return
end

% compute acceptable window limits
minTime = -abs(buffer);
maxTime = sessionDuration + abs(buffer);

% gather all times that must fall within the window
spikeTimes = sp.spike_times(:);
if ~isempty(spikeTimes) && (min(spikeTimes) < minTime || max(spikeTimes) > maxTime)
    error('glm:InvalidSpikeTimes', 'Spike times fall outside supported session bounds.');
end

if isempty(ev)
    return
end

allTimes = [ [ev.t_on], [ev.t_off] ];
if ~isempty(allTimes) && (min(allTimes) < minTime || max(allTimes) > maxTime)
    error('glm:InvalidLabels', 'Event times fall outside supported session bounds.');
end
end
