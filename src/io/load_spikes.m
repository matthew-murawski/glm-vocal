function sp = load_spikes(path)
%LOAD_SPIKES Load spike times for a single neuron from disk.
%   sp = LOAD_SPIKES(path) loads spike data from a MAT file and validates it.

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

data = load(path);
requiredFields = {'spike_times', 'neuron_id', 'session_id'};
for idx = 1:numel(requiredFields)
    field = requiredFields{idx};
    if ~isfield(data, field)
        error('glm:InvalidSpikesStruct', 'Missing required field %s in spike file.', field);
    end
end

sp = struct();
sp.spike_times = data.spike_times;
sp.neuron_id = data.neuron_id;
sp.session_id = data.session_id;

if ~isa(sp.spike_times, 'double') || ~isvector(sp.spike_times)
    error('glm:InvalidSpikeTimes', 'Spike times must be a double vector.');
end

sp.spike_times = sp.spike_times(:);
if any(~isfinite(sp.spike_times))
    error('glm:InvalidSpikeTimes', 'Spike times must be finite values.');
end

if numel(sp.spike_times) > 1 && any(diff(sp.spike_times) < 0)
    [sp.spike_times, sortIdx] = sort(sp.spike_times, 'ascend'); %#ok<ASGLU>
    warning('glm:UnsortedSpikeTimes', 'Spike times were unsorted and have been sorted in ascending order.');
end

end
