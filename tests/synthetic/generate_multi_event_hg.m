function data = generate_multi_event_hg(varargin)
%GENERATE_MULTI_EVENT_HG Create synthetic session with multiple event types.
%   Creates a synthetic high gamma session with multiple events of different
%   types (addressed, overheard, produced), each with distinct kernel shapes.
%   Tests the model's ability to fit multiple kernels simultaneously.
%
%   Usage:
%       data = generate_multi_event_hg()  % use defaults
%       data = generate_multi_event_hg('save', true)  % save to file
%
%   Output:
%       data - struct with fields:
%              * hg_data: high gamma trace struct
%              * events: event struct array (all events)
%              * stim: stimulus timeline and streams
%              * cfg: configuration
%              * ground_truth: true kernel parameters for all types

% parse inputs
p = inputParser;
p.addParameter('save', false, @islogical);
p.addParameter('output_path', 'tests/data/synthetic_multi_event.mat', @ischar);
parse(p, varargin{:});

% section parameters
duration_s = 20.0;
dt = 0.01;
fs = 100;

% baseline parameters
baseline_power = 10.0;
noise_std = 0.5;

% event timing (non-overlapping)
addressed_times = [3.0, 9.0, 15.0];
overheard_times = [5.0, 12.0];
produced_times = [7.0, 17.0];

% kernel parameters
% addressed: fast positive response
addressed_peak = 4.0;
addressed_delay = 0.1;
addressed_width = 0.15;

% overheard: slow positive response
overheard_peak = 2.0;
overheard_delay = 0.3;
overheard_width = 0.4;

% produced: biphasic (pre-dip and post-peak)
produced_predip_mag = -1.0;
produced_predip_time = -0.1;
produced_predip_width = 0.05;
produced_peak = 3.0;
produced_peak_time = 0.2;
produced_peak_width = 0.2;

% section build time base
t = (0:dt:duration_s-dt)';
nT = length(t);

% section create events
events = [];
event_idx = 1;

% heard_addressed events
for i = 1:length(addressed_times)
    events(event_idx).kind = 'perceived';
    events(event_idx).t_on = addressed_times(i);
    events(event_idx).t_off = addressed_times(i) + 0.5;
    events(event_idx).label = 'addressed';
    events(event_idx).heard_type = 'heard_addressed';
    event_idx = event_idx + 1;
end

% heard_overheard events
for i = 1:length(overheard_times)
    events(event_idx).kind = 'perceived';
    events(event_idx).t_on = overheard_times(i);
    events(event_idx).t_off = overheard_times(i) + 0.5;
    events(event_idx).label = 'overheard';
    events(event_idx).heard_type = 'heard_overheard';
    event_idx = event_idx + 1;
end

% produced events
for i = 1:length(produced_times)
    events(event_idx).kind = 'produced';
    events(event_idx).t_on = produced_times(i);
    events(event_idx).t_off = produced_times(i) + 0.5;
    events(event_idx).label = 'spontaneous';
    events(event_idx).produced_context = 'produced_spontaneous';
    event_idx = event_idx + 1;
end

% section build stimulus streams
stim = struct();
stim.t = t;
stim.dt = dt;

streams = struct();
streams.heard_any = zeros(nT, 1);
streams.heard_addressed = zeros(nT, 1);
streams.heard_overheard = zeros(nT, 1);
streams.produced_spontaneous = zeros(nT, 1);
streams.produced_after_heard = zeros(nT, 1);
streams.produced_after_produced = zeros(nT, 1);

% mark heard_addressed
for i = 1:length(addressed_times)
    idx = find(t >= addressed_times(i) & t < addressed_times(i) + 0.5);
    streams.heard_any(idx) = 1;
    streams.heard_addressed(idx) = 1;
end

% mark heard_overheard
for i = 1:length(overheard_times)
    idx = find(t >= overheard_times(i) & t < overheard_times(i) + 0.5);
    streams.heard_any(idx) = 1;
    streams.heard_overheard(idx) = 1;
end

% mark produced
for i = 1:length(produced_times)
    idx = find(t >= produced_times(i) & t < produced_times(i) + 0.5);
    streams.produced_spontaneous(idx) = 1;
end

% specify which heard fields to use separately
streams.heard_fields = {'heard_addressed', 'heard_overheard'};

stim.streams = streams;

% section build states
states = struct();
states.convo = zeros(nT, 1);
states.spon = ones(nT, 1);
stim.states = states;

% section generate high gamma response
power = baseline_power * ones(nT, 1);

% build kernels
heard_window = [0, 2.0];
produced_window = [-2.0, 3.0];

% addressed kernel (fast positive)
addressed_kernel_times = (heard_window(1):dt:heard_window(2))';
addressed_kernel = addressed_peak * exp(-((addressed_kernel_times - addressed_delay).^2) / (2 * addressed_width^2));

% overheard kernel (slow positive)
overheard_kernel_times = (heard_window(1):dt:heard_window(2))';
overheard_kernel = overheard_peak * exp(-((overheard_kernel_times - overheard_delay).^2) / (2 * overheard_width^2));

% produced kernel (biphasic)
produced_kernel_times = (produced_window(1):dt:produced_window(2))';
produced_kernel = produced_predip_mag * exp(-((produced_kernel_times - produced_predip_time).^2) / (2 * produced_predip_width^2)) + ...
                  produced_peak * exp(-((produced_kernel_times - produced_peak_time).^2) / (2 * produced_peak_width^2));

% convolve streams with kernels
for i = 1:nT
    % addressed responses
    if streams.heard_addressed(i) > 0
        for k = 1:length(addressed_kernel)
            future_idx = i + k - 1;
            if future_idx <= nT
                power(future_idx) = power(future_idx) + addressed_kernel(k);
            end
        end
    end

    % overheard responses
    if streams.heard_overheard(i) > 0
        for k = 1:length(overheard_kernel)
            future_idx = i + k - 1;
            if future_idx <= nT
                power(future_idx) = power(future_idx) + overheard_kernel(k);
            end
        end
    end

    % produced responses (includes pre-dip)
    if streams.produced_spontaneous(i) > 0
        % symmetric window includes negative lags
        for k = 1:length(produced_kernel)
            target_idx = i + k - 1 + round(produced_window(1)/dt);
            if target_idx > 0 && target_idx <= nT
                power(target_idx) = power(target_idx) + produced_kernel(k);
            end
        end
    end
end

% add noise
rng(123);  % different seed from single event
noise = noise_std * randn(nT, 1);
power = power + noise;

% ensure non-negative
power = max(power, 0);

% section build hg_data struct
hg_data = struct();
hg_data.power = power;
hg_data.fs = fs;
hg_data.t = t;
hg_data.session_id = 'synthetic_multi_event';
hg_data.channel_id = 1;

% section build configuration
cfg = struct();
cfg.dt = dt;
cfg.heard_window_s = heard_window;
cfg.produced_window_s = produced_window;
cfg.history_window_s = [0, 0.5];
cfg.heard_basis = struct('kind', 'raised_cosine');
cfg.produced_basis = struct('kind', 'raised_cosine');
cfg.preprocess.hg_binning_method = 'mean';
cfg.model.lambda = 0.05;  % moderate regularization
cfg.model.link = 'identity';
cfg.optimization.max_iter = 200;
cfg.optimization.tol_fun = 1e-6;
cfg.optimization.tol_grad = 1e-6;
cfg.exclude_predictors = {};

% section ground truth
ground_truth = struct();
ground_truth.baseline = baseline_power;
ground_truth.noise_std = noise_std;
ground_truth.duration_s = duration_s;
ground_truth.dt = dt;

% addressed kernel info
ground_truth.addressed.kernel_times = addressed_kernel_times;
ground_truth.addressed.kernel_response = addressed_kernel;
ground_truth.addressed.peak = addressed_peak;
ground_truth.addressed.delay = addressed_delay;
ground_truth.addressed.width = addressed_width;
ground_truth.addressed.event_times = addressed_times;

% overheard kernel info
ground_truth.overheard.kernel_times = overheard_kernel_times;
ground_truth.overheard.kernel_response = overheard_kernel;
ground_truth.overheard.peak = overheard_peak;
ground_truth.overheard.delay = overheard_delay;
ground_truth.overheard.width = overheard_width;
ground_truth.overheard.event_times = overheard_times;

% produced kernel info
ground_truth.produced.kernel_times = produced_kernel_times;
ground_truth.produced.kernel_response = produced_kernel;
ground_truth.produced.predip_mag = produced_predip_mag;
ground_truth.produced.predip_time = produced_predip_time;
ground_truth.produced.peak = produced_peak;
ground_truth.produced.peak_time = produced_peak_time;
ground_truth.produced.event_times = produced_times;

% section assemble output
data = struct();
data.hg_data = hg_data;
data.events = events;
data.stim = stim;
data.cfg = cfg;
data.ground_truth = ground_truth;

% section save if requested
if p.Results.save
    output_path = p.Results.output_path;
    output_dir = fileparts(output_path);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save(output_path, '-struct', 'data');
    fprintf('Saved synthetic multi-event session to: %s\n', output_path);
end

end
