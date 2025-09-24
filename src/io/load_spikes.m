function sp = load_spikes(path)
%LOAD_SPIKES Load spike times for a single neuron from disk.
%   sp = LOAD_SPIKES(path) loads spike data from a MAT file and validates it.

% validate the supplied file path before touching the filesystem
if nargin < 1
    error('glm:InvalidInput', 'Path to spike MAT file is required.');
end
if isstring(path)
    path = char(path);
end
if ~ischar(path) || isempty(path)
    error('glm:InvalidInput', 'Path must be a non-empty character vector or string scalar.');
end
if exist(path, 'file') ~= 2
    error('glm:FileNotFound', 'Spike file not found: %s', path);
end

% load the file payload and make sure the expected fields exist
data = load(path);
requiredFields = {'spike_times', 'neuron_id', 'session_id'};
for idx = 1:numel(requiredFields)
    field = requiredFields{idx};
    if ~isfield(data, field)
        error('glm:InvalidSpikesStruct', 'Missing required field %s in spike file.', field);
    end
end

% copy fields into a clean struct for downstream consumers
sp = struct();
sp.spike_times = data.spike_times;
sp.neuron_id = data.neuron_id;
sp.session_id = data.session_id;

% enforce double column vector spike times with finite entries
if ~isa(sp.spike_times, 'double') || ~isvector(sp.spike_times)
    error('glm:InvalidSpikeTimes', 'Spike times must be a double vector.');
end
sp.spike_times = sp.spike_times(:);
if any(~isfinite(sp.spike_times))
    error('glm:InvalidSpikeTimes', 'Spike times must be finite values.');
end

% tidy ordering so downstream histograms behave as expected
if numel(sp.spike_times) > 1 && any(diff(sp.spike_times) < 0)
    sp.spike_times = sort(sp.spike_times, 'ascend');
    warning('glm:UnsortedSpikeTimes', 'Spike times were unsorted and have been sorted in ascending order.');
end

end
