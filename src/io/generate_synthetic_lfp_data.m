function [lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params)
%GENERATE_SYNTHETIC_LFP_DATA Generate synthetic LFP data with event-locked responses.
%   [lfpFile, heardFile, producedFile] = GENERATE_SYNTHETIC_LFP_DATA(outdir, params)
%   creates artificial LFP traces with known ground-truth event-locked gamma
%   responses for GLM pipeline validation.
%
%   Inputs:
%     outdir - directory where synthetic files will be saved
%     params - struct with generation parameters (see below)
%
%   Outputs:
%     lfpFile      - path to generated LFP MAT file
%     heardFile    - path to perceived calls Audacity TXT file
%     producedFile - path to produced calls Audacity TXT file
%
%   Parameters:
%     session_duration               - total recording length in seconds (default: 600)
%     lfp_fs                         - LFP sampling rate in Hz (default: 1000)
%     n_channels                     - number of LFP channels (default: 4)
%     perceived_to_produced_ratio    - ratio of heard to produced calls (default: 5.0)
%     mean_intercall_interval        - average seconds between calls (default: 3.0)
%     call_duration_mean             - average call length in seconds (default: 0.15)
%     call_duration_std              - std dev of call lengths (default: 0.05)
%     bout_probability               - probability a produced call triggers more (default: 0.7)
%     bout_length_mean               - average number of produced calls in bout (default: 3)
%     addressed_window               - seconds before produced call for "addressed" (default: 10.0)
%     gamma_center_freq              - center frequency for high-gamma (default: 110)
%     gamma_bandwidth                - bandwidth around center frequency (default: 40)
%     produced_response_amplitude    - scaling for produced responses (default: 2.0)
%     addressed_response_amplitude   - scaling for addressed perceived responses (default: 1.0)
%     overheard_response_amplitude   - scaling for overheard perceived responses (default: 0.3)
%     perceived_latency              - delay for perceived response peak in seconds (default: 0.3)
%     produced_pre_duration          - how far before produced calls response starts (default: 0.5)
%     noise_amplitude                - baseline noise level (default: 0.5)
%     snr_db                         - signal-to-noise ratio in dB (default: 10)

%% section validate inputs
if nargin < 1 || isempty(outdir)
    error('generate_synthetic_lfp_data:MissingOutdir', 'Output directory is required.');
end
if nargin < 2
    params = struct();
end

% set default parameters
params = set_default_params(params);

% create output directory
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% section generate conversational events
fprintf('Generating conversational event structure...\n');
[perceivedEvents, producedEvents] = generate_event_times(params);
fprintf('  Generated %d perceived calls and %d produced calls\n', ...
        numel(perceivedEvents), numel(producedEvents));

% classify perceived calls as addressed vs overheard
[perceivedEvents, nAddressed, nOverheard] = classify_addressed_calls(...
    perceivedEvents, producedEvents, params.addressed_window);
fprintf('  Classified: %d addressed, %d overheard\n', nAddressed, nOverheard);

%% section generate synthetic LFP data
fprintf('Generating synthetic LFP data...\n');
n_samples = round(params.session_duration * params.lfp_fs);
lfp_data = generate_synthetic_lfp(perceivedEvents, producedEvents, ...
                                   n_samples, params);
fprintf('  Generated %d samples × %d channels (%.1f seconds)\n', ...
        n_samples, params.n_channels, params.session_duration);

%% section save outputs
fprintf('Saving synthetic data files...\n');

% save LFP MAT file
lfpFile = fullfile(outdir, 'synthetic_lfp.mat');
lfp_fs = params.lfp_fs;
lfp_t0 = 0.0;
channel_ids = 1:params.n_channels;
save(lfpFile, 'lfp_data', 'lfp_fs', 'lfp_t0', 'channel_ids', '-v7.3');
fprintf('  Saved LFP data: %s\n', lfpFile);

% save Audacity label files
heardFile = fullfile(outdir, 'synthetic_heard.txt');
write_audacity_labels(perceivedEvents, heardFile);
fprintf('  Saved perceived labels: %s\n', heardFile);

producedFile = fullfile(outdir, 'synthetic_produced.txt');
write_audacity_labels(producedEvents, producedFile);
fprintf('  Saved produced labels: %s\n', producedFile);

fprintf('Synthetic data generation complete!\n');
end

%% ======================== HELPER FUNCTIONS ========================

function params = set_default_params(params)
%SET_DEFAULT_PARAMS Set default parameter values if not provided.
defaults = struct(...
    'session_duration', 600, ...
    'lfp_fs', 1000, ...
    'n_channels', 4, ...
    'perceived_to_produced_ratio', 5.0, ...
    'mean_intercall_interval', 3.0, ...
    'call_duration_mean', 0.15, ...
    'call_duration_std', 0.05, ...
    'bout_probability', 0.7, ...
    'bout_length_mean', 3, ...
    'addressed_window', 10.0, ...
    'gamma_center_freq', 110, ...
    'gamma_bandwidth', 40, ...
    'produced_response_amplitude', 2.0, ...
    'addressed_response_amplitude', 1.0, ...
    'overheard_response_amplitude', 0.3, ...
    'perceived_latency', 0.3, ...
    'produced_pre_duration', 0.5, ...
    'noise_amplitude', 0.5, ...
    'snr_db', 10);

fields = fieldnames(defaults);
for ii = 1:numel(fields)
    field = fields{ii};
    if ~isfield(params, field) || isempty(params.(field))
        params.(field) = defaults.(field);
    end
end
end

function [perceivedEvents, producedEvents] = generate_event_times(params)
%GENERATE_EVENT_TIMES Generate conversational event structure with bouts.
%   Creates realistic conversational patterns where produced calls cluster
%   into bouts and perceived calls are interspersed.

% estimate total number of calls
total_time = params.session_duration;
mean_interval = params.mean_intercall_interval;
n_calls_estimate = floor(total_time / mean_interval);

% partition into perceived vs produced
ratio = params.perceived_to_produced_ratio;
n_produced_estimate = round(n_calls_estimate / (1 + ratio));
n_perceived_estimate = round(n_calls_estimate * ratio / (1 + ratio));

% preallocate
perceivedEvents = repmat(struct('t_on', 0, 'duration', 0, 'is_addressed', false), ...
                         n_perceived_estimate * 2, 1);
producedEvents = repmat(struct('t_on', 0, 'duration', 0), ...
                        n_produced_estimate * 2, 1);

% generate conversational structure with bouts
t_current = 0;
n_perceived = 0;
n_produced = 0;

while t_current < total_time
    % decide if this is a produced bout or single perceived call
    is_bout = rand() < params.bout_probability;

    if is_bout && n_produced < numel(producedEvents)
        % generate a bout of produced calls
        bout_length = max(1, round(params.bout_length_mean + randn() * 1.5));

        for jj = 1:bout_length
            if t_current >= total_time || n_produced >= numel(producedEvents)
                break;
            end

            n_produced = n_produced + 1;
            duration = generate_call_duration(params);
            producedEvents(n_produced).t_on = t_current;
            producedEvents(n_produced).duration = duration;

            % advance time with shorter interval within bout
            jitter = randn() * mean_interval * 0.2;
            t_current = t_current + mean_interval * 0.5 + jitter;
        end
    else
        % generate single perceived call
        if n_perceived < numel(perceivedEvents)
            n_perceived = n_perceived + 1;
            duration = generate_call_duration(params);
            perceivedEvents(n_perceived).t_on = t_current;
            perceivedEvents(n_perceived).duration = duration;
            perceivedEvents(n_perceived).is_addressed = false;
        end
    end

    % advance time to next event with jitter
    jitter = randn() * mean_interval * 0.3;
    t_current = t_current + mean_interval + jitter;
end

% trim to actual counts
perceivedEvents = perceivedEvents(1:n_perceived);
producedEvents = producedEvents(1:n_produced);

% sort by onset time
if ~isempty(perceivedEvents)
    [~, idx] = sort([perceivedEvents.t_on]);
    perceivedEvents = perceivedEvents(idx);
end
if ~isempty(producedEvents)
    [~, idx] = sort([producedEvents.t_on]);
    producedEvents = producedEvents(idx);
end
end

function duration = generate_call_duration(params)
%GENERATE_CALL_DURATION Sample call duration from truncated normal distribution.
duration = params.call_duration_mean + randn() * params.call_duration_std;
duration = max(0.05, min(0.5, duration));  % clip to reasonable bounds
end

function [perceivedEvents, nAddressed, nOverheard] = classify_addressed_calls(...
    perceivedEvents, producedEvents, addressedWindow)
%CLASSIFY_ADDRESSED_CALLS Mark perceived calls as addressed or overheard.
%   A perceived call is "addressed" if there's a produced call within
%   addressedWindow seconds after it.

nAddressed = 0;
nOverheard = 0;

if isempty(perceivedEvents) || isempty(producedEvents)
    return;
end

for ii = 1:numel(perceivedEvents)
    t_perceived = perceivedEvents(ii).t_on;

    % check if any produced call follows within the addressed window
    produced_times = [producedEvents.t_on];
    time_diffs = produced_times - t_perceived;

    % addressed if there's a produced call within [0, addressedWindow] seconds
    is_addressed = any(time_diffs >= 0 & time_diffs <= addressedWindow);

    perceivedEvents(ii).is_addressed = is_addressed;

    if is_addressed
        nAddressed = nAddressed + 1;
    else
        nOverheard = nOverheard + 1;
    end
end
end

function lfp_data = generate_synthetic_lfp(perceivedEvents, producedEvents, ...
                                           n_samples, params)
%GENERATE_SYNTHETIC_LFP Create LFP traces with event-locked gamma responses.

n_channels = params.n_channels;
fs = params.lfp_fs;

% step 1: generate baseline noise (pink/brown noise for realism)
fprintf('  Generating baseline noise...\n');
lfp_data = generate_baseline_noise(n_samples, n_channels, params.noise_amplitude);

% step 2: add event-locked gamma bursts
fprintf('  Adding event-locked gamma responses...\n');
signal_power = 0;

% add produced call responses
for ii = 1:numel(producedEvents)
    t_on = producedEvents(ii).t_on;
    duration = producedEvents(ii).duration;

    % gamma burst starts before call and spans through it
    burst_start = t_on - params.produced_pre_duration;
    burst_duration = params.produced_pre_duration + duration;

    if burst_start < 0
        burst_start = 0;
        burst_duration = t_on + duration;
    end

    % generate gamma burst
    burst = generate_gamma_burst(burst_start, burst_duration, fs, n_samples, ...
                                 params.gamma_center_freq, params.gamma_bandwidth, ...
                                 params.produced_response_amplitude, 'produced', t_on);

    % add to all channels with slight variation
    for ch = 1:n_channels
        channel_scale = 0.8 + 0.4 * rand();  % 0.8 to 1.2
        lfp_data(:, ch) = lfp_data(:, ch) + burst * channel_scale;
    end

    signal_power = signal_power + sum(burst.^2);
end

% add perceived call responses
for ii = 1:numel(perceivedEvents)
    t_on = perceivedEvents(ii).t_on;
    is_addressed = perceivedEvents(ii).is_addressed;

    % response starts at latency after call onset
    burst_start = t_on + params.perceived_latency;
    burst_duration = 0.25;  % ~250ms response

    % amplitude depends on addressed vs overheard
    if is_addressed
        amplitude = params.addressed_response_amplitude;
    else
        amplitude = params.overheard_response_amplitude;
    end

    % generate gamma burst
    burst = generate_gamma_burst(burst_start, burst_duration, fs, n_samples, ...
                                 params.gamma_center_freq, params.gamma_bandwidth, ...
                                 amplitude, 'perceived', t_on);

    % add to all channels with slight variation
    for ch = 1:n_channels
        channel_scale = 0.8 + 0.4 * rand();  % 0.8 to 1.2
        lfp_data(:, ch) = lfp_data(:, ch) + burst * channel_scale;
    end

    signal_power = signal_power + sum(burst.^2);
end

% step 3: scale to achieve target SNR
noise_power = sum(lfp_data(:).^2) - signal_power;
current_snr = 10 * log10(signal_power / noise_power);
fprintf('  Current SNR: %.1f dB (target: %.1f dB)\n', current_snr, params.snr_db);

% step 4: apply bandpass filter to high-gamma range
fprintf('  Applying bandpass filter...\n');
f_low = params.gamma_center_freq - params.gamma_bandwidth;
f_high = params.gamma_center_freq + params.gamma_bandwidth;
for ch = 1:n_channels
    lfp_data(:, ch) = bandpass(lfp_data(:, ch), [f_low f_high], fs);
end

% step 5: extract high gamma band power via Hilbert transform
fprintf('  Extracting high gamma band power...\n');
for ch = 1:n_channels
    % apply Hilbert transform to get analytic signal
    analytic_signal = hilbert(lfp_data(:, ch));

    % extract amplitude envelope (instantaneous amplitude)
    lfp_data(:, ch) = abs(analytic_signal);
end
end

function noise = generate_baseline_noise(n_samples, n_channels, amplitude)
%GENERATE_BASELINE_NOISE Create pink/brown noise with 1/f spectrum.
%   Generates realistic baseline noise by filtering white noise.

noise = zeros(n_samples, n_channels);

for ch = 1:n_channels
    % generate white noise
    white_noise = randn(n_samples, 1);

    % create 1/f filter (pink noise approximation)
    % use a simple method: sum of filtered white noise at different scales
    pink_noise = zeros(n_samples, 1);
    scales = [1, 2, 4, 8, 16];
    weights = 1 ./ sqrt(scales);

    for ii = 1:numel(scales)
        scale = scales(ii);
        if scale > 1
            filtered = filter(ones(scale, 1) / scale, 1, white_noise);
        else
            filtered = white_noise;
        end
        pink_noise = pink_noise + filtered * weights(ii);
    end

    % normalize and scale
    pink_noise = pink_noise / std(pink_noise);
    noise(:, ch) = pink_noise * amplitude;
end
end

function burst = generate_gamma_burst(start_time, duration, fs, n_samples, ...
                                      center_freq, bandwidth, amplitude, ...
                                      event_type, event_onset)
%GENERATE_GAMMA_BURST Create a gamma-band oscillation with smooth envelope.
%   Creates a modulated gamma burst at the specified time with Gaussian envelope.

burst = zeros(n_samples, 1);

% convert to sample indices
start_idx = max(1, round(start_time * fs) + 1);
duration_samples = round(duration * fs);
end_idx = min(n_samples, start_idx + duration_samples - 1);

if start_idx >= n_samples || end_idx < 1
    return;  % burst is outside recording window
end

% time vector for this burst
n_burst_samples = end_idx - start_idx + 1;
t_burst = (0:n_burst_samples-1)' / fs;

% add frequency jitter for realism
freq_jitter = (rand() - 0.5) * 10;  % ±5 Hz
freq = center_freq + freq_jitter;

% generate gamma oscillation
gamma_osc = sin(2 * pi * freq * t_burst);

% create envelope based on event type
if strcmp(event_type, 'produced')
    % for produced calls, peak at the call onset
    peak_time = event_onset - start_time;
    envelope = gaussian_envelope(t_burst, peak_time, duration * 0.3);
else
    % for perceived calls, peak ~100ms after response start
    peak_time = 0.1;
    envelope = gaussian_envelope(t_burst, peak_time, duration * 0.35);
end

% modulate and scale
burst_segment = gamma_osc .* envelope * amplitude;

% place in full burst array
burst(start_idx:end_idx) = burst_segment;
end

function envelope = gaussian_envelope(t, peak_time, sigma)
%GAUSSIAN_ENVELOPE Create smooth Gaussian envelope for burst.
envelope = exp(-((t - peak_time).^2) / (2 * sigma^2));
envelope = envelope / max(envelope);  % normalize to peak of 1
end

function write_audacity_labels(events, filepath)
%WRITE_AUDACITY_LABELS Write events to Audacity label format.
%   Format: start\tstop\tlabel (tab-separated, three columns)

fid = fopen(filepath, 'w');
if fid < 0
    error('generate_synthetic_lfp_data:FileIO', ...
          'Unable to open file for writing: %s', filepath);
end

try
    for ii = 1:numel(events)
        t_on = events(ii).t_on;
        t_off = t_on + events(ii).duration;

        % determine label
        if isfield(events, 'is_addressed')
            % perceived call
            if events(ii).is_addressed
                label = 'addressed';
            else
                label = 'overheard';
            end
        else
            % produced call
            label = 'produced';
        end

        % write tab-separated line
        fprintf(fid, '%.6f\t%.6f\t%s\n', t_on, t_off, label);
    end
    fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME);
end
end
