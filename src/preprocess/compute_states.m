function states = compute_states(ev, stim, stateCfg)
% build conversational and spontaneous state flags from labelled events.
%
% policy (produce-anchored, two-step lookback/forward with rollback):
%   1) build produced (target) sequences by merging produced calls with gaps
%      <= max_seq_gap_s (default 3.0 s).
%   2) for each produced sequence with onset P_on:
%        - look back respWin (e.g., 5 s). if heard onsets exist, take the
%          latest within that window and walk backward through heard onsets
%          that fall within backChainWin (0.5 s) hops; the earliest link in
%          that chain marks the conversational window start.
%        - let P_off be this produced sequence's offset. look forward 5.0 s:
%          if there is NO heard onset in (P_off, P_off + response_window_s],
%          end the window at P_off and stop for this seed.
%          otherwise, take the FIRST heard onset H1_on after P_off within 5 s,
%          then look forward 5.0 s from that heard call's offset H1_off:
%            * if there is NO produced onset in (H1_off, H1_off + response_window_s],
%              end the window (retroactively) at P_off and stop.
%            * if there IS a produced onset P_next_on within 5 s, extend the
%              window through the corresponding produced sequence's offset
%              P_next_off, set current produced = that sequence, and repeat
%              the forward checks (chain continues while deadlines are met).
%   3) if multiple produced sequences generate windows that overlap or touch,
%      merge them into one interval before rasterizing.
%
% config (stateCfg):
%   - response_window_s : reply window in seconds (e.g., 5.0)
%   - max_seq_gap_s     : max allowed same-actor gap for produced sequences
%                         (optional; default 3.0)
%
% outputs:
%   states.convo : logical mask on stim grid for conversational state
%   states.spon  : logical mask for spontaneous (complement within mask.good)

%% validate inputs
if nargin < 3
    error('glm:InvalidInput', 'events, stimulus, and state config structs are required.');
end
if ~isstruct(stim) || ~isfield(stim, 't') || ~isfield(stim, 'dt')
    error('glm:InvalidInput', 'stimulus struct must contain fields t and dt.');
end

t  = stim.t(:);
dt = stim.dt;
if isempty(t)
    states = struct('convo', false(0,1), 'spon', false(0,1));
    return
end
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('glm:InvalidInput', 'stimulus dt must be a positive finite scalar.');
end
nBins     = numel(t);
gridStart = t(1);

%% thresholds
respWin = state_double(stateCfg, 'response_window_s');     % e.g., 5.0
if isfield(stateCfg, 'max_seq_gap_s')
    maxGap = state_double(stateCfg, 'max_seq_gap_s');      % e.g., 3.0
else
    maxGap = 3.0;
end
backChainWin = 0.5;  % backward heard-chain hop in seconds
tol = 1e-9;  % small time tolerance for boundary comparisons

%% extract and sort events
heard    = select_events(ev, 'perceived');
produced = select_events(ev, 'produced');

% heard calls as individual events
H_on  = double([heard.t_on]);  H_on  = H_on(:);
H_off = double([heard.t_off]); H_off = H_off(:);
if ~isempty(H_on)
    [H_on, hOrd] = sort(H_on);
    H_off = H_off(hOrd);
else
    H_off = zeros(0,1);
end

% build produced sequences
P_seq = build_produced_sequences(produced, maxGap);
P_on  = double([P_seq.start_on]); P_on  = P_on(:);
P_off = double([P_seq.end_off]);  P_off = P_off(:);

%% derive conversational intervals using the produced-anchored rules
intervals = zeros(0,2);
for p = 1:numel(P_on)
    Pon  = P_on(p);
    Poff = P_off(p);

    % step 1: look back respWin for heard onset(s) and chain backward with 0.5 s hops.
    hBackIdx = find(H_on >= (Pon - respWin - tol) & H_on <= (Pon + tol));
    if isempty(hBackIdx)
        continue
    end

    latestIdx = hBackIdx(end);
    chainIdx  = latestIdx;
    H0_on     = H_on(chainIdx);

    while true
        prevMask = H_on < (H_on(chainIdx) - tol) & ...
            H_on >= (H_on(chainIdx) - backChainWin - tol);
        if ~any(prevMask)
            break
        end
        chainIdx = find(prevMask, 1, 'last');
        H0_on    = H_on(chainIdx);
    end

    % initialize current produced pointer and end candidate
    currP = p;
    endAt = P_off(currP);

    % step 2 & 3: forward chaining with rollback rule
    while true
        % from current produced offset, look for first heard onset within respWin
        Pcurr_off = P_off(currP);
        hFwdIdx = find(H_on > Pcurr_off + tol & H_on <= Pcurr_off + respWin + tol, 1, 'first');
        if isempty(hFwdIdx)
            % no heard within window -> end at current produced offset
            endAt = max(endAt, Pcurr_off);
            break
        end

        % found first heard onset; take its offset
        H1_on  = H_on(hFwdIdx);
        H1_off = H_off(hFwdIdx);

        % after heard offset, look for a produced onset within respWin
        % (choose the first produced onset > H1_off)
        pFwdIdx = find(P_on > H1_off + tol & P_on <= H1_off + respWin + tol, 1, 'first');
        if isempty(pFwdIdx)
            % rollback: despite a heard onset after Pcurr_off, no produced onset
            % within 5 s of that heard offset -> end at the previous produced offset
            endAt = max(endAt, Pcurr_off);
            break
        end

        % extend through the corresponding produced sequence's offset and continue
        currP = pFwdIdx;
        endAt = max(endAt, P_off(currP));
        % loop continues: try to find next heard after this produced, etc.
    end

    intervals(end+1, :) = [H0_on, endAt]; %#ok<AGROW>
end

% merge overlapping/touching intervals
intervals = merge_intervals(intervals, tol);

%% rasterize conversational intervals to the stimulus grid
convoMask = mark_intervals(intervals, gridStart, dt, nBins);

%% apply valid-data mask and expose spontaneous as the complement
if isfield(stim, 'mask') && isstruct(stim.mask) && isfield(stim.mask, 'good')
    goodMask = logical(stim.mask.good(:));
    if numel(goodMask) ~= nBins
        error('glm:InvalidInput', 'stim.mask.good must align with stim.t.');
    end
    convoMask = convoMask & goodMask;
    sponMask  = goodMask & ~convoMask;
else
    sponMask = ~convoMask;
end

states = struct('convo', convoMask, 'spon', sponMask);
end

% -------------------------------------------------------------------------
% helpers
% -------------------------------------------------------------------------

function val = state_double(stateCfg, fieldName)
% coerce configuration scalar to a finite non-negative double value.
if ~isfield(stateCfg, fieldName)
    error('glm:InvalidInput', 'state config must include field %s.', fieldName);
end
raw = stateCfg.(fieldName);
if ~isscalar(raw) || ~isnumeric(raw)
    error('glm:InvalidInput', 'state config field %s must be a numeric scalar.', fieldName);
end
val = double(raw);
if ~isfinite(val) || val < 0
    error('glm:InvalidInput', 'state config field %s must be finite and non-negative.', fieldName);
end
end

function eventsOut = select_events(eventsIn, kind)
% filter events by kind, tolerating empty or missing inputs gracefully.
if isempty(eventsIn)
    eventsOut = repmat(struct('kind','', 't_on',[], 't_off',[]), 0, 1);
    return
end
if ~isstruct(eventsIn) || ~all(isfield(eventsIn, {'kind','t_on','t_off'}))
    error('glm:InvalidInput', 'event structs must include kind, t_on, and t_off fields.');
end

kindValues = cellfun(@char, {eventsIn.kind}, 'UniformOutput', false);
mask = strcmp(kindValues, kind);
if ~any(mask)
    template = eventsIn(1);
    template.kind  = '';
    template.t_on  = [];
    template.t_off = [];
    eventsOut = repmat(template, 0, 1);
else
    eventsOut = eventsIn(mask);
end
end

function P_seq = build_produced_sequences(produced, maxGap)
% merge produced calls into sequences when inter-call gap <= maxGap.
P_seq = repmat(struct('start_on',0, 'end_off',0), 0, 1);
if isempty(produced)
    return
end

% sort produced by onset
[~, ord] = sort(double([produced.t_on]));
prod = produced(ord);

% initialize first sequence
curr.start_on = double(prod(1).t_on);
curr.end_off  = double(prod(1).t_off);

for ii = 2:numel(prod)
    eOn  = double(prod(ii).t_on);
    eOff = double(prod(ii).t_off);
    gap  = eOn - curr.end_off;

    if gap <= maxGap
        % merge into current sequence
        if eOff > curr.end_off
            curr.end_off = eOff;
        end
    else
        % close current and start new
        if curr.end_off > curr.start_on
            P_seq(end+1,1) = curr; %#ok<AGROW>
        end
        curr.start_on = eOn;
        curr.end_off  = eOff;
    end
end

% append final sequence
if curr.end_off > curr.start_on
    P_seq(end+1,1) = curr;
end
end

function intervals = merge_intervals(intervals, tol)
% merge overlapping or touching intervals; treat |gap| <= tol as touching.
if isempty(intervals)
    return
end

% sort by start
[~, ord] = sort(intervals(:,1));
iv = intervals(ord, :);

merged = iv(1, :);
for k = 2:size(iv,1)
    prev = merged(end, :);
    cur  = iv(k, :);
    if cur(1) <= prev(2) + tol  % overlap or touch
        merged(end, 2) = max(prev(2), cur(2));
    else
        merged(end+1, :) = cur; %#ok<AGROW>
    end
end
intervals = merged;
end

function mask = mark_intervals(intervals, gridStart, dt, nBins)
% convert continuous intervals to a logical mask on the stimulus grid.

mask = false(nBins, 1);
if isempty(intervals)
    return
end

for ii = 1:size(intervals, 1)
    tOn  = intervals(ii, 1);
    tOff = intervals(ii, 2);

    if ~(isfinite(tOn) && isfinite(tOff)) || tOff <= tOn
        continue
    end

    % inclusive start, exclusive end mapping
    startIdx = floor((tOn  - gridStart) / dt) + 1;
    edgeTol  = max([eps(dt), eps(tOff), eps(tOff - gridStart)]);
    endEdge  = (tOff - gridStart) - edgeTol;
    endIdx   = ceil(endEdge / dt);

    if endIdx < startIdx
        continue
    end

    startIdx = max(startIdx, 1);
    endIdx   = min(endIdx, nBins);
    if startIdx > nBins || endIdx < 1
        continue
    end

    mask(startIdx:endIdx) = true;
end
end