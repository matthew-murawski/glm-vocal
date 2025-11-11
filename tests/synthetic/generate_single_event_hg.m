function data = generate_single_event_hg(varargin)
%GENERATE_SINGLE_EVENT_HG Create robust synthetic session with single perceived call.
%   Creates a synthetic high gamma session with one perceived call and a
%   known Gaussian kernel response embedded in baseline with added noise.
%   This is the primary test case for end-to-end pipeline validation.
%
%   Usage:
%       data = generate_single_event_hg()  % use defaults
%       data = generate_single_event_hg('save', true)  % save to file
%
%   Output:
%       data - struct with fields:
%              * hg_data: high gamma trace struct (for load_high_gamma format)
%              * events: event struct array
%              * stim: stimulus timeline and streams
%              * cfg: configuration
%              * ground_truth: true kernel parameters for validation

% parse inputs
p = inputParser;
p.addParameter('save', false, @islogical);
p.addParameter('output_path', 'tests/data/synthetic_single_event.mat', @ischar);
parse(p, varargin{:});

% section parameters (robust for step 6)
duration_s = 10.0;     % 10 seconds
dt = 0.01;             % 10 ms bins, 100 Hz
fs = 100;              % high gamma sampling rate

% baseline parameters
baseline_power = 10.0;  % higher baseline

% event parameters
event_time = 5.0;      % single perceived call at middle of session

% kernel parameters (Gaussian bump)
kernel_peak = 5.0;      % peak amplitude (larger for better recovery)
kernel_width = 0.2;     % standard deviation (seconds)
kernel_delay = 0.15;    % delay to peak (seconds)

% noise parameters
noise_std = 0.5;        % Gaussian noise standard deviation

% section build time base
t = (0:dt:duration_s-dt)';
nT = length(t);

% section create event
events = struct();
events.kind = 'perceived';
events.t_on = event_time;
events.t_off = event_time + 0.5;  % 500ms duration
events.label = 'addressed';
events.heard_type = 'heard_addressed';

% section build stimulus streams
% initialize all streams to zero
stim = struct();
stim.t = t;
stim.dt = dt;

streams = struct();
streams.heard_any = zeros(nT, 1);
streams.heard_addressed = zeros(nT, 1);
streams.heard_overheard = zeros(nT, 1);

% mark the perceived event
event_idx = find(t >= events.t_on & t < events.t_off);
streams.heard_any(event_idx) = 1;
streams.heard_addressed(event_idx) = 1;

% add minimal produced streams (empty)
streams.produced_spontaneous = zeros(nT, 1);
streams.produced_after_heard = zeros(nT, 1);
streams.produced_after_produced = zeros(nT, 1);

% store streams
stim.streams = streams;

% section build states
states = struct();
states.convo = zeros(nT, 1);  % all spontaneous for simplicity
states.spon = ones(nT, 1);
stim.states = states;

% section generate high gamma response
power = baseline_power * ones(nT, 1);

% build gaussian kernel response
% causal window: [0, 2s] for heard events
heard_window = [0, 2.0];
kernel_times = (heard_window(1):dt:heard_window(2))';
kernel_response = kernel_peak * exp(-((kernel_times - kernel_delay).^2) / (2 * kernel_width^2));

% convolve event stream with kernel to get response
for i = 1:nT
    if streams.heard_addressed(i) > 0
        % add kernel response starting at this time
        for k = 1:length(kernel_response)
            future_idx = i + k - 1;
            if future_idx <= nT
                power(future_idx) = power(future_idx) + kernel_response(k);
            end
        end
    end
end

% add gaussian noise
rng(42);  % set seed for reproducibility
noise = noise_std * randn(nT, 1);
power = power + noise;

% ensure non-negative
power = max(power, 0);

% section build hg_data struct (mimics load_high_gamma output)
hg_data = struct();
hg_data.power = power;
hg_data.fs = fs;
hg_data.t = t;
hg_data.session_id = 'synthetic_single_event';
hg_data.channel_id = 1;

% section build configuration
cfg = struct();

% time discretization
cfg.dt = dt;

% kernel windows
cfg.heard_window_s = heard_window;
cfg.produced_window_s = [-2.0, 3.0];  % not used in this test
cfg.history_window_s = [0, 0.5];      % not used for high gamma

% basis configuration (use default raised cosine)
cfg.heard_basis = struct('kind', 'raised_cosine');
cfg.produced_basis = struct('kind', 'raised_cosine');

% preprocessing
cfg.preprocess = struct();
cfg.preprocess.hg_binning_method = 'mean';

% model settings
cfg.model = struct();
cfg.model.lambda = 0.01;  % small regularization for stability
cfg.model.link = 'identity';

% optimization
cfg.optimization = struct();
cfg.optimization.max_iter = 200;
cfg.optimization.tol_fun = 1e-6;
cfg.optimization.tol_grad = 1e-6;

% exclusions (spike history already excluded in wrapper)
cfg.exclude_predictors = {};

% section ground truth
ground_truth = struct();
ground_truth.baseline = baseline_power;
ground_truth.kernel_peak = kernel_peak;
ground_truth.kernel_width = kernel_width;
ground_truth.kernel_delay = kernel_delay;
ground_truth.kernel_times = kernel_times;
ground_truth.kernel_response = kernel_response;
ground_truth.event_time = event_time;
ground_truth.event_type = 'heard_addressed';
ground_truth.noise_std = noise_std;
ground_truth.duration_s = duration_s;
ground_truth.dt = dt;

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

    % create directory if it doesn't exist
    output_dir = fileparts(output_path);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % save
    save(output_path, '-struct', 'data');
    fprintf('Saved synthetic session to: %s\n', output_path);
end

end
