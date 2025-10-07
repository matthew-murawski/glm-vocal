function buildSyntheticData(varargin)
% builds synthetic spike times and vocal events for glm validation.
%
% usage: buildSyntheticData('Duration', 600, 'RespondsToProduced', false)

% prepare the simulation parameters and parse name-value arguments so the generator is reproducible and configurable. this also anchors the default output directory relative to the repository root.
parser = inputParser;
parser.FunctionName = 'buildSyntheticData';
currentDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(currentDir);
repoRoot = fileparts(srcDir);
defaultOutputDir = fullfile(repoRoot, 'data', 'Label Files', 'SyntheticSession');
addParameter(parser, 'RespondsToHeard', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'RespondsToProduced', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'Duration', 900, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'positive'}));
addParameter(parser, 'OutputDir', defaultOutputDir, @(x) (ischar(x) && ~isempty(x)) || (isstring(x) && isscalar(x)));
parse(parser, varargin{:});
opts = parser.Results;
if isstring(opts.OutputDir)
    opts.OutputDir = char(opts.OutputDir);
end

rng(1234);
dt = 0.001;
baseline_rate_hz = 5;
call_duration_range_s = [0.2, 0.5];

if ~exist(opts.OutputDir, 'dir')
    mkdir(opts.OutputDir);
end

fprintf('generating event times...\n');
% generate the synthetic vocal events by drawing inter-call intervals until the session duration is filled. this keeps each call non-overlapping while letting the mix of heard and produced tokens emerge stochastically.
events = struct('onset', {}, 'offset', {}, 'type', {});
lastOffset = 0;
while lastOffset < opts.Duration
    interval_s = max(1.0, 4.5 + 1.5 * randn());
    onset_s = lastOffset + interval_s;
    if onset_s >= opts.Duration
        break;
    end
    evtDuration = call_duration_range_s(1) + diff(call_duration_range_s) * rand();
    offset_s = min(onset_s + evtDuration, opts.Duration);
    if rand() < 0.2
        evtType = 'produced';
    else
        evtType = 'heard';
    end
    events(end + 1) = struct('onset', onset_s, 'offset', offset_s, 'type', evtType); %#ok<AGROW>
    lastOffset = offset_s;
end

% define the heard and produced response kernels that modulate the firing rate in downstream steps. the delays and normalizations here mirror the qualitative expectations for the glm test cases.
heard_t = 0:dt:1;
heard_delay_s = 0.3;
heard_tau_s = 0.05;
heard_kernel = zeros(size(heard_t));
heard_mask = heard_t >= heard_delay_s;
heard_shifted = heard_t(heard_mask) - heard_delay_s;
heard_kernel(heard_mask) = heard_shifted .* exp(-heard_shifted / heard_tau_s);
if max(heard_kernel) > 0
    heard_kernel = heard_kernel / max(heard_kernel) * 25;
end

produced_t = -2.5:dt:1;
produced_mu_s = -1.5;
produced_sigma_s = 0.4;
produced_kernel = exp(-((produced_t - produced_mu_s) .^ 2) / (2 * produced_sigma_s ^ 2));
if max(produced_kernel) > 0
    produced_kernel = produced_kernel / max(produced_kernel) * 40;
end

% accumulate the baseline firing rate and add event-driven kernels with slight variability. this loop respects the response flags so callers can silence either influence.
t_session = 0:dt:opts.Duration;
rate_hz = baseline_rate_hz * ones(size(t_session));
for idx = 1:numel(events)
    evt = events(idx);
    if strcmp(evt.type, 'heard') && opts.RespondsToHeard
        abs_t = evt.onset + heard_t;
        valid = abs_t <= opts.Duration;
        kernel_vals = heard_kernel(valid) * (1 + 0.2 * randn());
        sample_idx = max(1, min(numel(t_session), floor(abs_t(valid) / dt) + 1));
        rate_hz(sample_idx) = rate_hz(sample_idx) + kernel_vals;
    elseif strcmp(evt.type, 'produced') && opts.RespondsToProduced
        abs_t = evt.onset + produced_t;
        valid = abs_t >= 0 & abs_t <= opts.Duration;
        kernel_vals = produced_kernel(valid) * (1 + 0.2 * randn());
        sample_idx = max(1, min(numel(t_session), floor(abs_t(valid) / dt) + 1));
        rate_hz(sample_idx) = rate_hz(sample_idx) + kernel_vals;
    end
end
rate_hz = max(rate_hz, 0);

fprintf('generating spikes...\n');
% sample spike counts from a poisson process and expand them back into a spike-time list. representing the spikes as timestamps keeps the output compatible with existing loaders.
expected_spikes_per_bin = rate_hz * dt;
spike_counts = poissrnd(expected_spikes_per_bin);
total_spikes = sum(spike_counts);
if total_spikes > 0
    spike_times = zeros(total_spikes, 1);
    nonzero_bins = find(spike_counts > 0);
    cursor = 1;
    for ii = 1:numel(nonzero_bins)
        bin_index = nonzero_bins(ii);
        count = spike_counts(bin_index);
        span = cursor:(cursor + count - 1);
        spike_times(span) = t_session(bin_index);
        cursor = cursor + count;
    end
else
    spike_times = zeros(0, 1);
end
sp.spike_times = spike_times;

fprintf('writing output files...\n');
% write spikes and event tables to disk so the existing pipeline can load them. both text files adopt the simple two-column layout used by the label importer.
spikeFile = fullfile(opts.OutputDir, 'synthetic_spike_times.mat');
neuron_id = 1; %#ok<NASGU>
session_id = 1; %#ok<NASGU>
spike_times = sp.spike_times; %#ok<NASGU>
save(spikeFile, 'neuron_id', 'session_id', 'spike_times');

if isempty(events)
    heard_events = struct('onset', {}, 'offset', {}, 'type', {});
    produced_events = heard_events;
else
    evt_types = {events.type};
    heard_events = events(strcmp(evt_types, 'heard'));
    produced_events = events(strcmp(evt_types, 'produced'));
end

writeEventFile(fullfile(opts.OutputDir, 'synthetic_heard.txt'), heard_events);
writeEventFile(fullfile(opts.OutputDir, 'synthetic_produced.txt'), produced_events);

fprintf('synthetic data written to %s\n', opts.OutputDir);
end

function writeEventFile(path, evtStruct)
% helper to serialize event lists into tab-delimited text files. wrapping this logic avoids repeating fopen and fclose calls in the main routine.
fid = fopen(path, 'w');
if fid == -1
    error('buildSyntheticData:FileOpenFailed', 'failed to open %s for writing', path);
end
cleanupObj = onCleanup(@() fclose(fid));
for ii = 1:numel(evtStruct)
    fprintf(fid, '%.6f\t%.6f\n', evtStruct(ii).onset, evtStruct(ii).offset);
end
end
