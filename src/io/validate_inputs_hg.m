function validation = validate_inputs_hg(hg)
%VALIDATE_INPUTS_HG Validate high gamma data prior to processing.
%   validation = VALIDATE_INPUTS_HG(hg) ensures that high gamma power
%   trace is numeric, finite, non-negative, and has variability. Time
%   vector must be monotonically increasing.
%
%   Returns a validation struct with fields:
%       valid       - true if all checks pass
%       warnings    - cell array of warning messages
%       errors      - cell array of error messages

% initialize validation struct
validation = struct();
validation.valid = true;
validation.warnings = {};
validation.errors = {};

% check that hg struct has required fields
requiredFields = {'power', 'fs', 't'};
for idx = 1:numel(requiredFields)
    field = requiredFields{idx};
    if ~isfield(hg, field)
        validation.errors{end+1} = sprintf('Missing required field: %s', field);
        validation.valid = false;
    end
end

% if missing critical fields, stop here
if ~validation.valid
    return
end

% validate power trace
power = hg.power;

% check type and shape
if ~isnumeric(power)
    validation.errors{end+1} = 'power must be numeric';
    validation.valid = false;
    return
end

if ~isvector(power)
    validation.errors{end+1} = 'power must be a vector';
    validation.valid = false;
    return
end

power = power(:);  % ensure column vector
nSamples = numel(power);

% check for finite values
if any(~isfinite(power))
    validation.errors{end+1} = 'power contains non-finite values (NaN or Inf)';
    validation.valid = false;
end

% check for negative values (warning, not error, since they can be clipped)
if any(power < 0)
    n_negative = sum(power < 0);
    validation.warnings{end+1} = sprintf(...
        'power contains %d negative values (%.2f%%). These should be set to zero.', ...
        n_negative, 100 * n_negative / nSamples);
end

% check for all zeros
if all(power == 0)
    validation.errors{end+1} = 'power is all zeros - no signal to model';
    validation.valid = false;
end

% check for variability
if std(power) < eps
    validation.errors{end+1} = sprintf(...
        'power has no variability (std = %g < eps). Cannot model constant data.', ...
        std(power));
    validation.valid = false;
end

% validate sampling rate
fs = hg.fs;
if ~isscalar(fs) || ~isfinite(fs)
    validation.errors{end+1} = 'fs must be a finite scalar';
    validation.valid = false;
elseif fs <= 0
    validation.errors{end+1} = sprintf('fs must be positive (got %g)', fs);
    validation.valid = false;
end

% validate time vector
t = hg.t;

% check type and shape
if ~isnumeric(t) || ~isvector(t)
    validation.errors{end+1} = 't must be a numeric vector';
    validation.valid = false;
    return
end

t = t(:);  % ensure column vector

% check length matches power
if numel(t) ~= numel(power)
    validation.errors{end+1} = sprintf(...
        't length (%d) does not match power length (%d)', ...
        numel(t), numel(power));
    validation.valid = false;
end

% check for finite values
if any(~isfinite(t))
    validation.errors{end+1} = 't contains non-finite values';
    validation.valid = false;
end

% check monotonically increasing
if numel(t) > 1
    dt = diff(t);
    if any(dt <= 0)
        validation.errors{end+1} = 't must be strictly monotonically increasing';
        validation.valid = false;
    end
end

end
