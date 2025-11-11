function hg = load_high_gamma(path, varargin)
%LOAD_HIGH_GAMMA Load high gamma power traces from MAT file.
%   hg = LOAD_HIGH_GAMMA(path) loads high gamma data from a MAT file
%   exported by export_continuous_hg_for_glm and validates it.
%
%   hg = LOAD_HIGH_GAMMA(path, 'channel_id', id) selects a single channel.
%
%   hg = LOAD_HIGH_GAMMA(path, 'reduce', 'mean') reduces across channels
%   using the specified method ('mean' or 'median').
%
%   Output struct fields:
%       power       - [nSamples × 1] double, non-negative power trace
%       fs          - scalar sampling rate (Hz)
%       t           - [nSamples × 1] time vector (seconds)
%       session_id  - string session identifier (if available)
%       channel_id  - selected channel ID or 'reduced'

% parse input arguments
p = inputParser;
addRequired(p, 'path', @(x) ischar(x) || isstring(x));
addParameter(p, 'channel_id', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'reduce', '', @(x) ischar(x) || isstring(x));
parse(p, path, varargin{:});

% validate the supplied file path before touching the filesystem
if isstring(path)
    path = char(path);
end
if ~ischar(path) || isempty(path)
    error('glm:InvalidInput', 'Path must be a non-empty character vector or string scalar.');
end
if exist(path, 'file') ~= 2
    error('glm:FileNotFound', 'High gamma file not found: %s', path);
end

% load the file payload and make sure the expected fields exist
data = load(path);
if ~isfield(data, 'lfp_data')
    error('glm:InvalidHGStruct', 'Missing required field lfp_data in high gamma file.');
end

lfp_data = data.lfp_data;
requiredFields = {'lfp_data', 'lfp_fs', 'lfp_t0', 'channel_ids'};
for idx = 1:numel(requiredFields)
    field = requiredFields{idx};
    if ~isfield(lfp_data, field)
        error('glm:InvalidHGStruct', 'Missing required field lfp_data.%s in high gamma file.', field);
    end
end

% extract raw power matrix and metadata
power_matrix = lfp_data.lfp_data;  % [nSamples × nChannels]
fs = lfp_data.lfp_fs;
t0 = lfp_data.lfp_t0;
channel_ids = lfp_data.channel_ids;

% validate dimensions and types
if ~isa(power_matrix, 'double') || ndims(power_matrix) > 2
    error('glm:InvalidHGData', 'lfp_data.lfp_data must be a 2D double matrix.');
end
if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
    error('glm:InvalidHGData', 'lfp_data.lfp_fs must be a positive scalar.');
end
if ~isscalar(t0) || ~isfinite(t0)
    error('glm:InvalidHGData', 'lfp_data.lfp_t0 must be a finite scalar.');
end

[nSamples, nChannels] = size(power_matrix);
if ~isvector(channel_ids) || numel(channel_ids) ~= nChannels
    error('glm:InvalidHGData', 'lfp_data.channel_ids must be a vector matching number of channels.');
end

% select single channel or reduce across channels
channel_id_arg = p.Results.channel_id;
reduce_method = p.Results.reduce;

if ~isempty(channel_id_arg) && ~isempty(reduce_method)
    error('glm:InvalidInput', 'Cannot specify both channel_id and reduce. Choose one.');
end

if ~isempty(channel_id_arg)
    % single channel selection
    chan_idx = find(channel_ids == channel_id_arg, 1);
    if isempty(chan_idx)
        error('glm:InvalidInput', 'Channel ID %d not found in data.', channel_id_arg);
    end
    power = power_matrix(:, chan_idx);
    selected_channel_id = channel_id_arg;
elseif ~isempty(reduce_method)
    % reduce across channels
    reduce_method = char(reduce_method);
    switch lower(reduce_method)
        case 'mean'
            power = mean(power_matrix, 2);
        case 'median'
            power = median(power_matrix, 2);
        otherwise
            error('glm:InvalidInput', 'Reduce method must be ''mean'' or ''median''.');
    end
    selected_channel_id = 'reduced';
else
    % default: use first channel
    if nChannels > 1
        warning('glm:MultipleChannels', ...
            'Multiple channels found. Using first channel (ID=%d). Specify channel_id or reduce to choose.', ...
            channel_ids(1));
    end
    power = power_matrix(:, 1);
    selected_channel_id = channel_ids(1);
end

% ensure column vector
power = power(:);

% handle negative values
if any(power < 0)
    n_negative = sum(power < 0);
    warning('glm:NegativeValues', ...
        'Found %d negative power values (%.2f%%). Setting to zero.', ...
        n_negative, 100 * n_negative / nSamples);
    power(power < 0) = 0;
end

% construct time vector: prefer metadata.timestamps if available
if isfield(data, 'metadata') && isfield(data.metadata, 'timestamps')
    t = data.metadata.timestamps;
    if ~isvector(t) || numel(t) ~= nSamples
        warning('glm:InvalidTimestamps', ...
            'metadata.timestamps length mismatch. Constructing from lfp_t0 and lfp_fs.');
        t = t0 + (0:nSamples-1)' / fs;
    else
        t = t(:);  % ensure column vector
    end
else
    % construct from t0 and fs
    t = t0 + (0:nSamples-1)' / fs;
end

% extract session_id if available
if isfield(data, 'metadata') && isfield(data.metadata, 'session_id')
    session_id = data.metadata.session_id;
else
    session_id = '';
end

% build output struct
hg = struct();
hg.power = power;
hg.fs = fs;
hg.t = t;
hg.session_id = session_id;
hg.channel_id = selected_channel_id;

end
