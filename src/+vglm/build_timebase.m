function t = vglm.build_timebase(session_dur, dt)
% build uniform time grid for glm bins.

% check scalar inputs for basic validity.
if ~(isscalar(session_dur) && isfinite(session_dur) && session_dur > 0)
    error('vglm:build_timebase:BadSessionDur', ...
        'session_dur must be a positive finite scalar.');
end
if ~(isscalar(dt) && isfinite(dt) && dt > 0)
    error('vglm:build_timebase:BadDt', ...
        'dt must be a positive finite scalar.');
end

% ensure the ratio is (nearly) integer so the grid tiles perfectly.
ratio = session_dur / dt;
if abs(ratio - round(ratio)) > 1e-9
    error('vglm:build_timebase:NonIntegerRatio', ...
        'session_dur/dt must be an integer within tolerance.');
end
nbins = round(ratio);

% build the column vector of bin start times.
starts = (0:(nbins - 1)) * dt;
t = starts(:);
end
