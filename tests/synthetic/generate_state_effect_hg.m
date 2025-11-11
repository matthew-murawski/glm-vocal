function data = generate_state_effect_hg(varargin)
%GENERATE_STATE_EFFECT_HG Create synthetic session with state-dependent baselines.
%   Creates a session with alternating conversational and spontaneous states,
%   where baseline power differs by state. Tests the model's ability to
%   estimate both kernel responses and state-dependent baseline shifts.
%
%   Usage:
%       data = generate_state_effect_hg()  % use defaults
%       data = generate_state_effect_hg('save', true)  % save to file
%
%   Output:
%       data - struct with fields:
%              * hg_data: high gamma trace struct
%              * events: event struct array
%              * stim: stimulus timeline and streams with states
%              * cfg: configuration
%              * ground_truth: true kernel and state coefficients

% parse inputs
p = inputParser;
p.addParameter('save', false, @islogical);
p.addParameter('output_path', 'tests/data/synthetic_state_effect.mat', @ischar);
parse(p, varargin{:});

% section parameters
duration_s = 30.0;
dt = 0.01;
fs = 100;

% baseline parameters by state
baseline_spon = 8.0;
baseline_convo = 12.0;  % +4 boost
baseline_boost_convo = baseline_convo - baseline_spon;  % +4

% state epochs
% [0-10s]: spontaneous
% [10-20s]: conversational
% [20-30s]: spontaneous
epoch_bounds = [0, 10, 20, 30];

% noise
noise_std = 0.5;

% event timing: 2 perceived calls per epoch
event_times = [3.0, 7.0,      % spontaneous epoch 1
               13.0, 17.0,    % conversational epoch
               23.0, 27.0];   % spontaneous epoch 2

% kernel parameters (same across all epochs)
kernel_peak = 3.0;
kernel_delay = 0.15;
kernel_width = 0.2;

% section build time base
t = (0:dt:duration_s-dt)';
nT = length(t);

% section create events
events = [];
for i = 1:length(event_times)
    events(i).kind = 'perceived';
    events(i).t_on = event_times(i);
    events(i).t_off = event_times(i) + 0.5;
    events(i).label = 'addressed';
    events(i).heard_type = 'heard_addressed';
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

% mark events
for i = 1:length(event_times)
    idx = find(t >= event_times(i) & t < event_times(i) + 0.5);
    streams.heard_any(idx) = 1;
    streams.heard_addressed(idx) = 1;
end

streams.heard_fields = {'heard_addressed'};
stim.streams = streams;

% section build states
states = struct();
states.convo = zeros(nT, 1);
states.spon = zeros(nT, 1);

% conversational state: [10-20s]
convo_idx = (t >= epoch_bounds(2)) & (t < epoch_bounds(3));
states.convo(convo_idx) = 1;

% spontaneous state: [0-10s] and [20-30s]
spon_idx = ((t >= epoch_bounds(1)) & (t < epoch_bounds(2))) | ...
           ((t >= epoch_bounds(3)) & (t <= epoch_bounds(4)));
states.spon(spon_idx) = 1;

stim.states = states;

% section generate high gamma response
power = zeros(nT, 1);

% add state-dependent baseline
power(spon_idx) = baseline_spon;
power(convo_idx) = baseline_convo;

% build kernel
heard_window = [0, 2.0];
kernel_times = (heard_window(1):dt:heard_window(2))';
kernel_response = kernel_peak * exp(-((kernel_times - kernel_delay).^2) / (2 * kernel_width^2));

% convolve events with kernel
for i = 1:nT
    if streams.heard_addressed(i) > 0
        for k = 1:length(kernel_response)
            future_idx = i + k - 1;
            if future_idx <= nT
                power(future_idx) = power(future_idx) + kernel_response(k);
            end
        end
    end
end

% add noise
rng(456);  % unique seed
noise = noise_std * randn(nT, 1);
power = power + noise;

% ensure non-negative
power = max(power, 0);

% section build hg_data struct
hg_data = struct();
hg_data.power = power;
hg_data.fs = fs;
hg_data.t = t;
hg_data.session_id = 'synthetic_state_effect';
hg_data.channel_id = 1;

% section build configuration
cfg = struct();
cfg.dt = dt;
cfg.heard_window_s = heard_window;
cfg.produced_window_s = [-2.0, 3.0];
cfg.history_window_s = [0, 0.5];
cfg.heard_basis = struct('kind', 'raised_cosine');
cfg.produced_basis = struct('kind', 'raised_cosine');
cfg.preprocess.hg_binning_method = 'mean';
cfg.model.lambda = 0.05;
cfg.model.link = 'identity';
cfg.optimization.max_iter = 200;
cfg.optimization.tol_fun = 1e-6;
cfg.optimization.tol_grad = 1e-6;
cfg.exclude_predictors = {};

% section ground truth
ground_truth = struct();
ground_truth.baseline_spon = baseline_spon;
ground_truth.baseline_convo = baseline_convo;
ground_truth.state_coeff_convo = baseline_boost_convo;  % +4
ground_truth.state_coeff_spon = 0;  % reference
ground_truth.noise_std = noise_std;
ground_truth.duration_s = duration_s;
ground_truth.dt = dt;
ground_truth.epoch_bounds = epoch_bounds;

% kernel info
ground_truth.kernel.kernel_times = kernel_times;
ground_truth.kernel.kernel_response = kernel_response;
ground_truth.kernel.peak = kernel_peak;
ground_truth.kernel.delay = kernel_delay;
ground_truth.kernel.width = kernel_width;
ground_truth.kernel.event_times = event_times;

% state timing
ground_truth.convo_mask = convo_idx;
ground_truth.spon_mask = spon_idx;

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
    fprintf('Saved synthetic state effect session to: %s\n', output_path);
end

end
