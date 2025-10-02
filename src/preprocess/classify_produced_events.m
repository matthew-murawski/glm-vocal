function categoryIdx = classify_produced_events(events, producedMask, lookbackWindowS, cfg)
% section purpose
% classify produced events either by conversational context or by acoustic call type.

if nargin < 4
    cfg = struct();
end
if nargin < 3 || isempty(lookbackWindowS)
    lookbackWindowS = 5.0;
end
if nargin < 2 || isempty(producedMask)
    producedMask = false(size(events));
end

mode = resolve_split_mode(cfg);

switch mode
    case 'context'
        categoryIdx = classify_by_context(events, producedMask, lookbackWindowS);
    case 'call_type'
        categoryIdx = classify_by_call_type(events, producedMask);
    otherwise
        error('classify_produced_events:InvalidMode', 'Unsupported produced split mode: %s', mode);
end

categoryIdx.mode = mode;
end

function categoryIdx = classify_by_context(events, producedMask, lookbackWindowS)
baseNames = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
categoryIdx = initialize_category_struct(baseNames);

if isempty(events)
    return
end

producedMask = logical(producedMask(:));
if ~any(producedMask)
    return
end

events = events(:);
tOnAll = arrayfun(@(ev) double(ev.t_on), events);
tOnAll = tOnAll(:);

kindCells = cellfun(@char, {events.kind}, 'UniformOutput', false);
kindCells = kindCells(:);
isRecognised = strcmp(kindCells, 'perceived') | strcmp(kindCells, 'produced');

producedIdx = find(producedMask);
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

function categoryIdx = classify_by_call_type(events, producedMask)
baseTypes = {'phee', 'twitter', 'trill', 'trillphee'};
producedMask = logical(producedMask(:));

if isempty(events) || ~any(producedMask)
    names = make_produced_fieldnames(baseTypes);
    categoryIdx = initialize_category_struct(names);
    return
end

events = events(:);
producedIdx = find(producedMask);
labels = cell(numel(producedIdx), 1);
for ii = 1:numel(producedIdx)
    evLabel = '';
    if isfield(events, 'label')
        evLabel = events(producedIdx(ii)).label;
    end
    labels{ii} = sanitize_call_type(evLabel);
end

observedTypes = unique(labels, 'stable');
allTypes = unique([baseTypes(:); observedTypes(:)], 'stable');
fieldNames = make_produced_fieldnames(allTypes);

categoryIdx = initialize_category_struct(fieldNames);
for kk = 1:numel(allTypes)
    matches = strcmp(labels, allTypes{kk});
    categoryIdx.(fieldNames{kk}) = producedIdx(matches);
end
end

function categoryIdx = initialize_category_struct(names)
categoryIdx = struct();
categoryIdx.names = names(:)';
for ii = 1:numel(names)
    categoryIdx.(names{ii}) = zeros(0, 1);
end
end

function mode = resolve_split_mode(cfg)
mode = 'context';
if ~isstruct(cfg)
    return
end
if ~isfield(cfg, 'produced_split_mode') || isempty(cfg.produced_split_mode)
    return
end

raw = lower(strtrim(string(cfg.produced_split_mode)));
if numel(raw) ~= 1
    error('classify_produced_events:InvalidMode', 'produced_split_mode must be a scalar string.');
end
candidate = char(raw);
if ~(strcmp(candidate, 'context') || strcmp(candidate, 'call_type'))
    error('classify_produced_events:InvalidMode', 'produced_split_mode must be ''context'' or ''call_type''.');
end
mode = candidate;
end

function fieldNames = make_produced_fieldnames(typeNames)
typeNames = cellfun(@(x) char(x), typeNames, 'UniformOutput', false);
fieldNames = cell(size(typeNames));
for ii = 1:numel(typeNames)
    fieldNames{ii} = ['produced_', sanitize_field_suffix(typeNames{ii})];
end
end

function clean = sanitize_field_suffix(name)
clean = lower(strtrim(char(name)));
clean = regexprep(clean, '[^a-z0-9]', '');
if isempty(clean)
    clean = 'unknown';
end
end

function out = sanitize_call_type(label)
if isa(label, 'string')
    label = char(label);
elseif ~ischar(label)
    label = '';
end
raw = lower(strtrim(label));
if isempty(raw)
    out = 'unknown';
    return
end
canonical = regexprep(raw, '[_\-\s]', '');
if contains(canonical, 'trillphee')
    out = 'trillphee';
elseif contains(canonical, 'trill')
    out = 'trill';
elseif contains(canonical, 'twitter')
    out = 'twitter';
elseif contains(canonical, 'phee')
    out = 'phee';
else
    tmp = regexprep(canonical, '\d+$', '');
    if isempty(tmp)
        tmp = canonical;
    end
    if isempty(tmp)
        tmp = 'unknown';
    end
    out = tmp;
end
end
