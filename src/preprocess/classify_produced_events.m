function categoryIdx = classify_produced_events(events, producedMask, lookbackWindowS)
% section purpose
% classify each produced event by looking back over recent heard/produced activity within a fixed window.

if nargin < 3 || isempty(lookbackWindowS)
    lookbackWindowS = 5.0;
end
if nargin < 2 || isempty(producedMask)
    producedMask = false(size(events));
end

categoryIdx = struct('produced_spontaneous', zeros(0, 1), ...
    'produced_after_heard', zeros(0, 1), ...
    'produced_after_produced', zeros(0, 1));

if isempty(events) || ~any(producedMask)
    return
end

producedMask = logical(producedMask(:));

events = events(:);
tOnAll = arrayfun(@(ev) double(ev.t_on), events);
tOnAll = tOnAll(:);

kindCells = cellfun(@char, {events.kind}, 'UniformOutput', false);
kindCells = kindCells(:);
isRecognised = strcmp(kindCells, 'perceived') | strcmp(kindCells, 'produced');

producedIdx = find(producedMask);
if isempty(producedIdx)
    return
end

spontaneous = false(size(producedIdx));
afterHeard = false(size(producedIdx));
afterProduced = false(size(producedIdx));

for ii = 1:numel(producedIdx)
    currIdx = producedIdx(ii);
    currOn = tOnAll(currIdx);
    if ~isfinite(currOn)
        spontaneous(ii) = true;
        continue
    end

    windowStart = currOn - lookbackWindowS;
    windowMask = isRecognised & isfinite(tOnAll) & tOnAll < currOn & tOnAll >= windowStart;
    if ~any(windowMask)
        spontaneous(ii) = true;
        continue
    end

    candidates = find(windowMask);
    tCandidates = tOnAll(candidates);
    [~, bestPos] = max(tCandidates(:));
    if isempty(bestPos) || ~isfinite(tCandidates(bestPos))
        spontaneous(ii) = true;
        continue
    end
    prevIdx = candidates(bestPos);

    prevKind = kindCells{prevIdx};
    if strcmp(prevKind, 'perceived')
        afterHeard(ii) = true;
    elseif strcmp(prevKind, 'produced')
        afterProduced(ii) = true;
    else
        spontaneous(ii) = true;
    end
end

fallback = ~(afterHeard | afterProduced | spontaneous);
spontaneous = spontaneous | fallback;

categoryIdx.produced_spontaneous = producedIdx(spontaneous);
categoryIdx.produced_after_heard = producedIdx(afterHeard);
categoryIdx.produced_after_produced = producedIdx(afterProduced);
end
