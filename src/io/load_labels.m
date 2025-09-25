function events = load_labels(path, defaultKind)
% load_labels load labelled events from mat or audacity txt files.
%
% events = load_labels(path) reads either a mat struct array (fields kind,
% onset, offset, label) or an audacity txt file formatted as
% "start\tstop\toptional_label". the output struct array normalizes fields
% to kind (char), t_on (double), t_off (double), and label (string). when
% parsing txt files, an optional default kind can be supplied to override
% missing or non-produced/perceived labels.

%% validate inputs
% we guard against missing paths and resolve supported string types.
if nargin < 1
    error('glm:InvalidInput', 'Path to labels file is required.');
end
path = ensure_char_path(path);
if exist(path, 'file') ~= 2
    error('glm:FileNotFound', 'Label file not found: %s', path);
end

%% sanitize optional default kind
% when provided, default kinds must be either produced or perceived.
if nargin < 2 || strlength(string(defaultKind)) == 0
    defaultKind = '';
else
    defaultKind = char(sanitize_kind(defaultKind));
end

%% dispatch based on file extension
% we currently support mat structs and audacity txt label exports.
[~, ~, ext] = fileparts(lower(path));
switch ext
    case '.mat'
        events = load_from_mat(path);
    case '.txt'
        events = load_from_txt(path, defaultKind);
    otherwise
        error('glm:UnsupportedLabelFormat', 'Labels must be provided as a MAT struct or an Audacity TXT file.');
end

%% normalize empty outputs to zero-length struct arrays
% downstream code expects an empty but typed struct when no events exist.
if isempty(events)
    events = repmat(empty_event_record(), 0, 1);
end
end

function events = load_from_txt(path, defaultKind)
% parse an audacity-style txt file into normalized event records.
fid = fopen(path, 'r');
if fid < 0
    error('glm:FileIO', 'Unable to open label file: %s', path);
end
cleaner = onCleanup(@() fclose(fid));

records = repmat(empty_event_record(), 0, 1);
lineIdx = 0;
while true
    rawLine = fgetl(fid);
    if ~ischar(rawLine)
        break
    end
    lineIdx = lineIdx + 1;
    line = strtrim(rawLine);
    if isempty(line) || startsWith(line, '#')
        continue
    end

    parts = strsplit(line, '\t');
    if numel(parts) < 2
        error('glm:InvalidLabelsStruct', 'TXT labels must contain onset and offset columns (line %d).', lineIdx);
    end

    tOnRaw = str2double(parts{1});
    tOffRaw = str2double(parts{2});
    if isnan(tOnRaw) || isnan(tOffRaw)
        error('glm:InvalidLabelTimes', 'TXT labels must use numeric onset/offset values (line %d).', lineIdx);
    end
    [tOn, tOff] = sanitize_times(tOnRaw, tOffRaw);

    if numel(parts) >= 3
        labelRaw = strtrim(strjoin(parts(3:end), '\t'));
    else
        labelRaw = '';
    end
    labelString = string(labelRaw);
    if strlength(labelString) == 0
        labelString = "";
    else
        labelString = labelString(1);
    end

    kindCandidate = '';
    if strlength(labelString) > 0
        try
            kindCandidate = sanitize_kind(labelString);
        catch %#ok<CTCH>
            kindCandidate = '';
        end
    end
    if isempty(kindCandidate)
        if isempty(defaultKind)
            error('glm:InvalidLabelKind', 'TXT labels require produced/perceived labels or a default kind (line %d).', lineIdx);
        end
        kindCandidate = defaultKind;
    end

    rec = empty_event_record();
    rec.kind = kindCandidate;
    rec.t_on = tOn;
    rec.t_off = tOff;
    rec.label = labelString;
    records(end+1, 1) = rec; %#ok<AGROW>
end

events = records;
end

function events = load_from_mat(path)
% pull the first struct variable from the mat file payload
raw = load(path);
structVars = fieldnames(raw);
selected = [];
for idx = 1:numel(structVars)
    candidate = raw.(structVars{idx});
    if isstruct(candidate)
        selected = candidate;
        break
    end
end
if isempty(selected) || ~isstruct(selected)
    error('glm:InvalidLabelsStruct', 'MAT file must contain a struct array named or containing event fields.');
end

% ensure the required event fields are present on every element
requiredFields = {'kind', 'onset', 'offset', 'label'};
for idx = 1:numel(requiredFields)
    if ~all(isfield(selected, requiredFields{idx}))
        error('glm:InvalidLabelsStruct', 'Event struct array must contain field: %s.', requiredFields{idx});
    end
end

% normalize each entry into the canonical event record shape
n = numel(selected);
events = repmat(empty_event_record(), n, 1);
for ii = 1:n
    entry = selected(ii);
    events(ii) = normalize_struct_record(entry);
end
end

function rec = normalize_struct_record(entry)
% convert the user supplied fields into normalized representations
kindStr = sanitize_kind(entry.kind);
[tOn, tOff] = sanitize_times(entry.onset, entry.offset);
labelStr = sanitize_label(entry.label);

% build a tidy record that downstream code can rely on
rec = empty_event_record();
rec.kind = kindStr;
rec.t_on = tOn;
rec.t_off = tOff;
rec.label = labelStr;
end

function kindStr = sanitize_kind(kindRaw)
% unwrap cells and standardize casing for the kind string
if iscell(kindRaw)
    if isempty(kindRaw)
        kindRaw = '';
    else
        kindRaw = kindRaw{1};
    end
end
kindStr = lower(strtrim(string(kindRaw)));

% validate the semantic meaning of the kind value
if strlength(kindStr) == 0
    error('glm:InvalidLabelKind', 'Event kind must be a non-empty string.');
end
if ~(kindStr == "produced" || kindStr == "perceived")
    error('glm:InvalidLabelKind', 'Event kind must be ''produced'' or ''perceived''.');
end
kindStr = char(kindStr);
end

function [tOn, tOff] = sanitize_times(tOnRaw, tOffRaw)
% check that provided times are scalars and numeric
if ~isscalar(tOnRaw) || ~isscalar(tOffRaw)
    error('glm:InvalidLabelTimes', 'Event times must be scalar values.');
end
if ~isnumeric(tOnRaw) || ~isnumeric(tOffRaw)
    error('glm:InvalidLabelTimes', 'Event times must be numeric.');
end

% guard against non-finite values and inverted intervals
if ~isfinite(tOnRaw) || ~isfinite(tOffRaw)
    error('glm:InvalidLabelTimes', 'Event times must be finite.');
end
tOn = double(tOnRaw);
tOff = double(tOffRaw);
if tOff < tOn
    error('glm:InvalidLabelTimes', 'Event end time must be greater than or equal to start time.');
end
end

function labelStr = sanitize_label(labelRaw)
% unwrap cells and coerce labels to string scalars
if iscell(labelRaw)
    if isempty(labelRaw)
        labelRaw = "";
    else
        labelRaw = labelRaw{1};
    end
end
labelStr = string(labelRaw);

% trim whitespace and default to empty when nothing provided
if isempty(labelStr)
    labelStr = "";
else
    labelStr = strtrim(labelStr(1));
end
end

function rec = empty_event_record()
% provide a consistent empty record template
rec = struct('kind', '', 't_on', 0, 't_off', 0, 'label', "");
end

function pathOut = ensure_char_path(pathIn)
% coerce supported string types into a char vector
if isstring(pathIn)
    if numel(pathIn) ~= 1
        error('glm:InvalidInput', 'Path must be a single string scalar.');
    end
    pathOut = char(pathIn);
elseif ischar(pathIn)
    pathOut = pathIn;
else
    error('glm:InvalidInput', 'Path must be a character vector or string scalar.');
end
end
