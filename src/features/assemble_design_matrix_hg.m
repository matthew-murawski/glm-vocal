function [X, y, colmap] = assemble_design_matrix_hg(hg_binned, events, stim, cfg)
%ASSEMBLE_DESIGN_MATRIX_HG Build design matrix for high gamma GLM.
%   [X, y, colmap] = ASSEMBLE_DESIGN_MATRIX_HG(hg_binned, events, stim, cfg)
%   creates the design matrix for modeling continuous high gamma power.
%   This is a thin wrapper around assemble_design_matrix that excludes
%   spike history regressors (not applicable for continuous power).
%
%   Inputs:
%       hg_binned - struct from bin_high_gamma with fields:
%                   * power: [n × 1] binned power trace
%                   * t: [n × 1] time vector matching stim.t
%       events    - struct array with processed events (produced/perceived)
%       stim      - struct with time base and event streams
%       cfg       - configuration with kernel settings
%
%   Outputs:
%       X      - [n × p] sparse design matrix
%       y      - [n × 1] power trace (response variable)
%       colmap - struct mapping column indices to regressor names

% section input validation
if ~isstruct(hg_binned) || ~isfield(hg_binned, 'power')
    error('glm:InvalidInput', 'hg_binned must be a struct with field: power');
end

if ~isstruct(stim) || ~isfield(stim, 't')
    error('glm:InvalidInput', 'stim must be a struct with field: t');
end

if ~isstruct(cfg)
    error('glm:InvalidInput', 'cfg must be a struct');
end

% section extract response variable
y = hg_binned.power(:);
n = length(y);

% validate response variable
if ~isnumeric(y) || ~isreal(y)
    error('glm:InvalidInput', 'power must be numeric and real');
end

if any(~isfinite(y))
    error('glm:InvalidInput', 'power contains NaN or Inf values');
end

if length(stim.t) ~= n
    error('glm:DimensionMismatch', ...
        'stim.t length (%d) must match binned power length (%d)', ...
        length(stim.t), n);
end

% section modify config to exclude spike history
cfg_hg = cfg;

% ensure exclude_predictors field exists
if ~isfield(cfg_hg, 'exclude_predictors')
    cfg_hg.exclude_predictors = {};
elseif ischar(cfg_hg.exclude_predictors)
    cfg_hg.exclude_predictors = {cfg_hg.exclude_predictors};
end

% add spike_history to exclusions
if ~ismember('spike_history', cfg_hg.exclude_predictors)
    cfg_hg.exclude_predictors{end+1} = 'spike_history';
end

% override spike history settings for safety
if isfield(cfg_hg, 'features')
    cfg_hg.features.include_spike_history = false;
    cfg_hg.features.spike_history_window = [];
end

% section call existing infrastructure
% pass empty spike counts (third argument) since we're using continuous power
empty_spikes = zeros(n, 1);

% check if events is a struct or if we need to extract streams
if isstruct(stim) && isfield(stim, 'streams')
    streams = stim.streams;
else
    error('glm:InvalidInput', 'stim must contain streams field with event streams');
end

% check if states are provided
if isstruct(stim) && isfield(stim, 'states')
    states = stim.states;
else
    % create default states if not provided
    states = struct();
    states.convo = zeros(n, 1);
    states.spon = ones(n, 1);
end

% call the existing design matrix assembler
% signature: assemble_design_matrix(streams, states, sps, cfg, stim)
Xd = assemble_design_matrix(streams, states, empty_spikes, cfg_hg, stim);

X = Xd.X;
colmap = Xd.colmap;

% section validate consistency
if size(X, 1) ~= n
    error('glm:DimensionMismatch', ...
        'Design matrix rows (%d) do not match response length (%d)', ...
        size(X, 1), n);
end

% check for NaN or Inf in design matrix
if any(~isfinite(X(:)))
    error('glm:InvalidOutput', 'Design matrix contains NaN or Inf values');
end

% verify no spike history columns present
if isfield(colmap, 'spike_history')
    error('glm:InvalidOutput', 'Design matrix incorrectly contains spike_history columns');
end

% section validate design matrix structure
validate_design_matrix_hg(X, y, colmap);

end


function validate_design_matrix_hg(X, y, colmap)
%VALIDATE_DESIGN_MATRIX_HG Validate design matrix structure for high gamma.
%   Performs comprehensive checks on the design matrix, response variable,
%   and column mapping to ensure consistency and expected structure.

% check dimensions match
if size(X, 1) ~= length(y)
    error('glm:ValidationFailed', ...
        'Design matrix rows (%d) must match response length (%d)', ...
        size(X, 1), length(y));
end

% verify no spike history columns
if isfield(colmap, 'spike_history')
    error('glm:ValidationFailed', ...
        'Design matrix should not contain spike_history columns for high gamma GLM');
end

% check sparsity is reasonable
sparsity_pct = 100 * (1 - nnz(X) / numel(X));
if sparsity_pct < 50
    warning('glm:LowSparsity', ...
        'Design matrix sparsity is only %.1f%%. Expected >50%% for typical GLM.', ...
        sparsity_pct);
end

% verify expected regressor types are present
% at minimum, we should have intercept
if ~isfield(colmap, 'intercept')
    warning('glm:MissingIntercept', ...
        'Design matrix does not contain intercept column');
end

% check for heard or produced regressors
has_heard = isfield(colmap, 'heard_fields') && ~isempty(colmap.heard_fields);
has_produced = isfield(colmap, 'produced_fields') && ~isempty(colmap.produced_fields);

if ~has_heard && ~has_produced
    warning('glm:NoEventRegressors', ...
        'Design matrix contains no heard or produced event regressors');
end

% verify column indices are within bounds
p = size(X, 2);
fields = fieldnames(colmap);
for i = 1:length(fields)
    field = fields{i};
    if isstruct(colmap.(field)) && isfield(colmap.(field), 'cols')
        cols = colmap.(field).cols;
        if any(cols < 1) || any(cols > p)
            error('glm:ValidationFailed', ...
                'Column mapping for %s contains out-of-bounds indices', field);
        end
    end
end

end
