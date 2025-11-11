function filepath = generate_simple_fit_problem(varargin)
%GENERATE_SIMPLE_FIT_PROBLEM Generate synthetic GLM fitting problem.
%   filepath = GENERATE_SIMPLE_FIT_PROBLEM() creates a simple synthetic
%   problem for testing GLM fitting with known ground truth.
%
%   Optional parameters:
%       'n'         - number of time bins (default: 100)
%       'p'         - number of features (default: 5)
%       'snr'       - signal-to-noise ratio (default: 3)
%       'filepath'  - output path (default: tests/data/synthetic_fit.mat)

% parse input arguments
parser = inputParser;
addParameter(parser, 'n', 100, @(x) isnumeric(x) && x > 0);
addParameter(parser, 'p', 5, @(x) isnumeric(x) && x > 0);
addParameter(parser, 'snr', 3, @(x) isnumeric(x) && x > 0);
addParameter(parser, 'filepath', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

n = parser.Results.n;
p = parser.Results.p;
snr = parser.Results.snr;
filepath = parser.Results.filepath;

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

    filepath = fullfile(data_dir, 'synthetic_fit.mat');
end

% set random seed for reproducibility
rng(12345);

% create design matrix
% include intercept and p-1 random features
X = [ones(n, 1), randn(n, p-1)];

% create true weights with varying magnitudes
w_true = zeros(p, 1);
w_true(1) = 10;  % intercept (baseline)
for i = 2:p
    w_true(i) = (-1)^i * (p - i + 1) * 0.5;  % alternating signs, decreasing magnitude
end

% generate signal
signal = X * w_true;

% determine noise level from SNR
signal_power = norm(signal)^2 / n;
noise_power = signal_power / (snr^2);
noise_std = sqrt(noise_power);

% generate noise
noise = noise_std * randn(n, 1);

% generate observations
y = signal + noise;

% ensure non-negative (high gamma power must be positive)
y(y < 0) = 0;

% create second-difference penalty matrix
D = sparse(p-2, p);
for i = 1:p-2
    D(i, i) = 1;
    D(i, i+1) = -2;
    D(i, i+2) = 1;
end

% store problem data
problem = struct();
problem.X = X;
problem.y = y;
problem.w_true = w_true;
problem.D = D;
problem.n = n;
problem.p = p;
problem.snr = snr;
problem.noise_std = noise_std;
problem.signal = signal;
problem.noise = noise;

% save to file
save(filepath, 'problem', '-v7.3');

fprintf('Generated simple fitting problem:\n');
fprintf('  Time bins (n):   %d\n', n);
fprintf('  Features (p):    %d\n', p);
fprintf('  SNR:             %.2f\n', snr);
fprintf('  Noise std:       %.4f\n', noise_std);
fprintf('  Signal range:    [%.2f, %.2f]\n', min(signal), max(signal));
fprintf('  Data range:      [%.2f, %.2f]\n', min(y), max(y));
fprintf('\n');
fprintf('  True weights:\n');
for i = 1:p
    fprintf('    w[%d] = %7.3f\n', i, w_true(i));
end
fprintf('\n');
fprintf('  Saved to: %s\n', filepath);

end
