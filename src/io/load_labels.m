function events = load_labels(path)
%LOAD_LABELS Load labelled events from a MAT file.
%   events = LOAD_LABELS(path) expects a MAT file containing a struct array
%   with fields: kind, onset, offset, label. The output is a struct array
%   with fields kind (char vector), t_on (double), t_off (double), and
%   label (string) after validation and normalization.

% check input value and resolve to a usable path
if nargin < 1
    error('glm:InvalidInput', 'Path to labels file is required.');
end
path = ensure_char_path(path);
if exist(path, 'file') ~= 2
    error('glm:FileNotFound', 'Label file not found: %s', path);
end

% enforce mat-only inputs before loading the file contents
[~, ~, ext] = fileparts(lower(path));
if ~strcmp(ext, '.mat')
    error('glm:UnsupportedLabelFormat', 'Labels must be provided as a MAT file containing a struct array.');
end

% decode the events and fall back to an empty array when nothing is stored
events = load_from_mat(path);
if isempty(events)
    events = repmat(empty_event_record(), 0, 1);
end
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
    pathOut = char(pathIn);
elseif ischar(pathIn)
    pathOut = pathIn;
else
    error('glm:InvalidInput', 'Path must be a character vector or string scalar.');
end
end
