function stim = build_timebase(ev, sp, dt)
%BUILD_TIMEBASE Placeholder constructor for the stimulus timebase.
%   stim = BUILD_TIMEBASE(ev, sp, dt) returns an empty scaffold.

unused = {ev, sp}; %#ok<NASGU>
stim = struct('t', [], 'dt', dt, 'mask', struct('good', []));
end
