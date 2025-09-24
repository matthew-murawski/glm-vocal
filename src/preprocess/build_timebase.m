function stim = build_timebase(ev, sp, dt)
%BUILD_TIMEBASE Construct an evenly spaced time axis for the session.
%   stim = BUILD_TIMEBASE(ev, sp, dt) infers the session duration from the
%   latest event or spike, then builds a column vector of bin centers and a
%   logical mask marking valid bins.

% validate the requested step size before using it
if nargin < 3 || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('glm:InvalidTimebase', 'Time step dt must be a positive finite scalar.');
end

% gather all relevant times from spikes and events
lastSpike = infer_last_time_from_spikes(sp);
lastEvent = infer_last_time_from_events(ev);

% choose the maximum observed time to cover the full session
T = max([0, lastSpike, lastEvent]);

% compute the number of bins needed to include the final observation
numBins = floor(T / dt + 0.5) + 1;

% generate the time vector and corresponding mask
stim = struct();
stim.t = (0:numBins-1)' * dt;
stim.dt = dt;
stim.mask = struct('good', true(numBins, 1));
end

function lastSpike = infer_last_time_from_spikes(sp)
% derive the latest spike time, guarding against empty inputs
if isempty(sp) || ~isstruct(sp) || ~isfield(sp, 'spike_times') || isempty(sp.spike_times)
    lastSpike = 0;
    return
end
spikeTimes = sp.spike_times(:);
if isempty(spikeTimes)
    lastSpike = 0;
else
    lastSpike = spikeTimes(end);
end
end

function lastEvent = infer_last_time_from_events(ev)
% obtain the furthest event offset if available
if isempty(ev) || ~isstruct(ev) || ~isfield(ev, 't_off') || isempty([ev.t_off])
    lastEvent = 0;
    return
end
lastEvent = max([ev.t_off]);
end
