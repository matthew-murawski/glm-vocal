function tools = struct_tools()
%STRUCT_TOOLS Utility accessors for struct operations.
%   tools = STRUCT_TOOLS() returns function handles for merge/get/set helpers.

% expose the public api as a bundle of function handles
tools = struct( ...
    'merge', @struct_merge, ...
    'get', @struct_get, ...
    'set', @struct_set);
end

function out = struct_merge(base, override)
%STRUCT_MERGE Deep merge of two structs, favoring override fields.
%   out = STRUCT_MERGE(base, override) recursively merges structs.

% normalise empty inputs so downstream logic is simpler
if nargin < 1 || isempty(base)
    base = struct();
end
if nargin < 2 || isempty(override)
    override = struct();
end

% fall back to plain replacement when either side is not a struct
if ~isstruct(base) || ~isstruct(override)
    out = override;
    return
end

% walk override fields and merge recursively when both sides are structs
out = base;
overrideFields = fieldnames(override);
for idx = 1:numel(overrideFields)
    field = overrideFields{idx};
    valueOverride = override.(field);
    if isfield(base, field)
        valueBase = base.(field);
    else
        valueBase = [];
    end

    if isstruct(valueBase) && isstruct(valueOverride)
        out.(field) = struct_merge(valueBase, valueOverride);
    else
        out.(field) = valueOverride;
    end
end
end

function v = struct_get(s, path, default)
%STRUCT_GET Retrieve a nested field using dot notation.
%   v = STRUCT_GET(s, path, default) returns default when any field is missing.

% handle missing path/default arguments up front
if nargin < 3
    default = [];
end
if nargin < 2 || isempty(path)
    v = s;
    return
end

% split the path and traverse until we fall off or succeed
segments = split_path(path);
current = s;
for idx = 1:numel(segments)
    key = segments{idx};
    if ~isstruct(current) || ~isfield(current, key)
        v = default;
        return
    end
    current = current.(key);
end
v = current;
end

function s = struct_set(s, path, value)
%STRUCT_SET Set a nested field using dot notation, creating structs as needed.
%   s = STRUCT_SET(s, path, value) materializes intermediate structs.

% default to empty struct or replace entirely when no path provided
if nargin < 1 || isempty(s)
    s = struct();
end
if nargin < 2 || isempty(path)
    s = value;
    return
end

% peel off the head segment and recurse through the struct tree
segments = split_path(path);
if numel(segments) == 1
    s.(segments{1}) = value;
    return
end

head = segments{1};
tailPath = strjoin(segments(2:end), '.');
if ~isfield(s, head) || ~isstruct(s.(head))
    s.(head) = struct();
end
s.(head) = struct_set(s.(head), tailPath, value);
end

function parts = split_path(path)
%SPLIT_PATH Split a dot-delimited path into components.

% coerce supported string types and validate input
if isstring(path)
    path = char(path);
end
if ~ischar(path)
    error('struct_tools:InvalidPath', 'Path must be a character vector or string scalar.');
end

% split into parts, preserving empty strings when someone passes ''
path = strtrim(path);
if isempty(path)
    parts = {''};
    return
end
parts = strsplit(path, '.');
end
