function lfp = load_lfp(path)
%LOAD_LFP Load continuous LFP data from disk.
%   lfp = LOAD_LFP(path) loads LFP data from a MAT file and validates it.
%
%   Expected fields in MAT file:
%     lfp_data     - (n_samples × n_channels) or (n_channels × n_samples) double matrix
%     lfp_fs       - scalar sampling rate in Hz
%     lfp_t0       - scalar start time in seconds (for alignment with spike times)
%     channel_ids  - (optional) vector or cell array of channel identifiers
%
%   Returns:
%     lfp.data        - (n_samples × n_channels) double matrix
%     lfp.fs          - sampling rate (Hz)
%     lfp.t0          - start time (seconds)
%     lfp.channel_ids - channel identifiers (1:n_channels if not provided)
%     lfp.n_samples   - number of time samples
%     lfp.n_channels  - number of channels
%     lfp.duration    - recording duration in seconds

% validate the supplied file path before touching the filesystem
if nargin < 1
    error('glm:InvalidInput', 'Path to LFP MAT file is required.');
end
if isstring(path)
    path = char(path);
end
if ~ischar(path) || isempty(path)
    error('glm:InvalidInput', 'Path must be a non-empty character vector or string scalar.');
end
if exist(path, 'file') ~= 2
    error('glm:FileNotFound', 'LFP file not found: %s', path);
end

% load the file payload and make sure the expected fields exist
data = load(path);
requiredFields = {'lfp_data', 'lfp_fs', 'lfp_t0'};
for idx = 1:numel(requiredFields)
    field = requiredFields{idx};
    if ~isfield(data, field)
        error('glm:InvalidLFPStruct', 'Missing required field %s in LFP file.', field);
    end
end

% copy fields into a clean struct for downstream consumers
lfp = struct();

% validate and normalize lfp_data to (n_samples × n_channels)
if ~isa(data.lfp_data, 'double') || ~ismatrix(data.lfp_data)
    error('glm:InvalidLFPData', 'LFP data must be a double matrix.');
end
if any(~isfinite(data.lfp_data(:)))
    error('glm:InvalidLFPData', 'LFP data contains non-finite values (NaN or Inf).');
end

% determine orientation and standardize to (n_samples × n_channels)
[dim1, dim2] = size(data.lfp_data);
if dim1 > dim2
    % assume (n_samples × n_channels) - typical case
    lfp.data = data.lfp_data;
    lfp.n_samples = dim1;
    lfp.n_channels = dim2;
else
    % assume (n_channels × n_samples) - transpose to standard form
    lfp.data = data.lfp_data';
    lfp.n_samples = dim2;
    lfp.n_channels = dim1;
    warning('glm:LFPTransposed', 'LFP data was transposed from (n_channels × n_samples) to (n_samples × n_channels).');
end

% validate sampling rate
if ~isscalar(data.lfp_fs) || ~isa(data.lfp_fs, 'double') || data.lfp_fs <= 0
    error('glm:InvalidSamplingRate', 'Sampling rate (lfp_fs) must be a positive scalar.');
end
lfp.fs = data.lfp_fs;

% validate start time
if ~isscalar(data.lfp_t0) || ~isa(data.lfp_t0, 'double') || ~isfinite(data.lfp_t0)
    error('glm:InvalidStartTime', 'Start time (lfp_t0) must be a finite scalar.');
end
lfp.t0 = data.lfp_t0;

% handle optional channel_ids
if isfield(data, 'channel_ids')
    if isnumeric(data.channel_ids)
        lfp.channel_ids = data.channel_ids(:);
    elseif iscell(data.channel_ids)
        lfp.channel_ids = data.channel_ids(:);
    else
        error('glm:InvalidChannelIDs', 'channel_ids must be a numeric vector or cell array.');
    end

    if numel(lfp.channel_ids) ~= lfp.n_channels
        error('glm:ChannelIDMismatch', 'Number of channel_ids (%d) does not match number of channels (%d).', ...
              numel(lfp.channel_ids), lfp.n_channels);
    end
else
    % default to 1:n_channels
    lfp.channel_ids = (1:lfp.n_channels)';
end

% compute derived quantities
lfp.duration = (lfp.n_samples - 1) / lfp.fs;

end
