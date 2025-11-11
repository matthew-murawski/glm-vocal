function filepath = generate_simple_hg_session(varargin)
%GENERATE_SIMPLE_HG_SESSION Generate minimal synthetic high gamma session.
%   filepath = GENERATE_SIMPLE_HG_SESSION() creates a flat (constant) high
%   gamma power trace with all required fields and saves to
%   tests/data/synthetic_simple.mat.
%
%   Optional parameters:
%       'duration'      - session duration in seconds (default: 10)
%       'fs'            - sampling rate in Hz (default: 100)
%       'baseline'      - constant power value (default: 10.0)
%       'n_channels'    - number of channels (default: 1)
%       'filepath'      - output path (default: tests/data/synthetic_simple.mat)
%       'dt'            - bin width for stim.t (default: 0.01 s)

% parse input arguments
p = inputParser;
addParameter(p, 'duration', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'baseline', 10.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'n_channels', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'filepath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'dt', 0.01, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

duration = p.Results.duration;
fs = p.Results.fs;
baseline = p.Results.baseline;
n_channels = p.Results.n_channels;
filepath = p.Results.filepath;
dt = p.Results.dt;

% determine output path
if isempty(filepath)
    % get the directory containing this script
    script_path = mfilename('fullpath');
    [script_dir, ~, ~] = fileparts(script_path);

    % construct path to tests/data/
    data_dir = fullfile(script_dir, '..', 'data');

    % ensure directory exists
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    filepath = fullfile(data_dir, 'synthetic_simple.mat');
end

% generate time vector and flat power trace
nSamples = round(duration * fs);
t0 = 0.0;

% create flat power matrix [nSamples × nChannels]
power_matrix = baseline * ones(nSamples, n_channels);

% create channel IDs
channel_ids = (1:n_channels)';

% construct metadata timestamps
timestamps = t0 + (0:nSamples-1)' / fs;

% build lfp_data struct (matches export_continuous_hg_for_glm format)
lfp_data = struct();
lfp_data.lfp_data = power_matrix;
lfp_data.lfp_fs = fs;
lfp_data.lfp_t0 = t0;
lfp_data.channel_ids = channel_ids;

% build metadata struct
metadata = struct();
metadata.timestamps = timestamps;
metadata.session_id = 'synthetic_simple';
metadata.duration = duration;
metadata.n_channels = n_channels;
metadata.baseline = baseline;
metadata.description = 'Flat constant power trace for testing I/O';

% build stim struct with time base for binning
nBins = floor(duration / dt);
stim = struct();
stim.t = (0:nBins-1)' * dt;
stim.dt = dt;

% save to file
save(filepath, 'lfp_data', 'metadata', 'stim', '-v7.3');

fprintf('Generated simple synthetic HG session:\n');
fprintf('  Duration: %.1f s\n', duration);
fprintf('  Sampling rate: %.1f Hz\n', fs);
fprintf('  Baseline: %.2f\n', baseline);
fprintf('  Channels: %d\n', n_channels);
fprintf('  Samples: %d\n', nSamples);
fprintf('  Bin width (dt): %.4f s\n', dt);
fprintf('  Bins: %d\n', nBins);
fprintf('  Saved to: %s\n', filepath);

end
