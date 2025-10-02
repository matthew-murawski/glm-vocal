function eventsOut = consolidate_twitter_bouts(eventsIn, boutWindow)
% consolidate_twitter_bouts collapse rapid twitter syllables into bouts.
%
% eventsOut = consolidate_twitter_bouts(eventsIn, boutWindow) returns a copy
% of the events array where produced-twitter calls occurring within
% boutWindow seconds of a bout onset are merged into a single event located
% at the first syllable onset. non-twitter events are passed through
% unchanged.

if nargin < 2 || isempty(boutWindow)
    boutWindow = 1.5;
end

if isempty(eventsIn)
    eventsOut = eventsIn;
    return
end

if ~isstruct(eventsIn) || ~all(isfield(eventsIn, {'kind', 't_on', 't_off'}))
    error('consolidate_twitter_bouts:InvalidEvents', 'event structs must include kind, t_on, and t_off fields.');
end

window = double(boutWindow);
if ~isscalar(window) || ~isfinite(window) || window <= 0
    error('consolidate_twitter_bouts:InvalidWindow', 'bout window must be a positive finite scalar.');
end

events = eventsIn(:);
kindCells = cellfun(@char, {events.kind}, 'UniformOutput', false);
kindCells = kindCells(:);
producedMask = strcmpi(kindCells, 'produced');
producedMask = producedMask(:);

if ~isfield(events, 'label')
twitterMask = false(numel(events), 1);
else
    labels = {events.label};
    twitterMask = false(numel(events), 1);
    for ii = 1:numel(labels)
        lbl = labels{ii};
        if isa(lbl, 'string')
            lbl = char(lbl);
        end
        if ~ischar(lbl)
            continue
        end
        canonical = lower(strtrim(lbl));
        canonical = regexprep(canonical, '[_\-\s]', '');
        twitterMask(ii) = contains(canonical, 'twitter');
    end
end

twitterMask = twitterMask(:);
targetMask = producedMask & twitterMask;
targetMask = targetMask(:);
if ~any(targetMask)
    eventsOut = eventsIn;
    return
end

twitterEvents = events(targetMask);
nonTwitterEvents = events(~targetMask);

[~, order] = sort(double([twitterEvents.t_on]));
twitterSorted = twitterEvents(order);

kept = repmat(empty_event_record(), 0, 1);
lastBoutStart = NaN;
for ii = 1:numel(twitterSorted)
    onset = double(twitterSorted(ii).t_on);
    if ~isfinite(onset)
        continue
    end
    if isnan(lastBoutStart) || (onset - lastBoutStart) > window
        lastBoutStart = onset;
        kept(end+1, 1) = twitterSorted(ii); %#ok<AGROW>
    end
end

combined = [nonTwitterEvents; kept];
[~, orderAll] = sort(double([combined.t_on]));
eventsOrdered = combined(orderAll);

% preserve original orientation (column)
eventsOrdered = eventsOrdered(:);

if size(eventsIn, 1) >= size(eventsIn, 2)
    eventsOut = eventsOrdered;
else
    eventsOut = reshape(eventsOrdered, 1, []);
end
end

function rec = empty_event_record()
rec = struct('kind', '', 't_on', [], 't_off', [], 'label', "");
end
