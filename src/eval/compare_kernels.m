function metrics = compare_kernels(fitted_kernels, true_kernels, kernel_names)
%COMPARE_KERNELS Compare fitted kernels to ground truth.
%   metrics = COMPARE_KERNELS(fitted_kernels, true_kernels, kernel_names)
%   computes comparison metrics for each kernel type and returns results
%   in a structured format with formatted table output.
%
%   Inputs:
%       fitted_kernels - struct with kernel fields from extract_kernels
%       true_kernels   - struct with ground truth kernel info
%       kernel_names   - cell array of kernel names to compare
%
%   Output:
%       metrics - struct with fields for each kernel:
%                 * correlation: correlation between fitted and true
%                 * peak_error_s: absolute error in peak timing (seconds)
%                 * peak_error_bins: peak error in bins
%                 * magnitude_ratio: fitted_peak / true_peak
%                 * magnitude_error: relative error in magnitude

% validate inputs
if ~isstruct(fitted_kernels)
    error('glm:InvalidInput', 'fitted_kernels must be a struct');
end

if ~isstruct(true_kernels)
    error('glm:InvalidInput', 'true_kernels must be a struct');
end

if ~iscell(kernel_names)
    kernel_names = {kernel_names};
end

% initialize output
metrics = struct();

% section compute metrics for each kernel
for i = 1:length(kernel_names)
    kernel_name = kernel_names{i};

    % check if kernel exists in both structs
    if ~isfield(fitted_kernels, kernel_name)
        warning('glm:MissingKernel', 'Fitted kernels missing: %s', kernel_name);
        continue;
    end

    if ~isfield(true_kernels, kernel_name)
        warning('glm:MissingKernel', 'True kernels missing: %s', kernel_name);
        continue;
    end

    fitted = fitted_kernels.(kernel_name);
    true_info = true_kernels.(kernel_name);

    % extract weights and times
    fitted_weights = fitted.weights(:);
    fitted_times = fitted.times(:);

    true_weights = true_info.kernel_response(:);
    true_times = true_info.kernel_times(:);

    % ensure same length (should be guaranteed by design)
    if length(fitted_weights) ~= length(true_weights)
        warning('glm:LengthMismatch', ...
            'Kernel %s: fitted length %d != true length %d', ...
            kernel_name, length(fitted_weights), length(true_weights));
        continue;
    end

    % correlation
    corr_val = corr(fitted_weights, true_weights);

    % find peaks
    [fitted_peak, fitted_peak_idx] = max(fitted_weights);
    [true_peak, true_peak_idx] = max(true_weights);

    fitted_peak_time = fitted_times(fitted_peak_idx);
    true_peak_time = true_times(true_peak_idx);

    % peak timing error
    peak_error_s = abs(fitted_peak_time - true_peak_time);

    % infer dt from time vector
    if length(fitted_times) > 1
        dt = median(diff(fitted_times));
    else
        dt = 0.01;  % default
    end
    peak_error_bins = peak_error_s / dt;

    % peak magnitude
    if true_peak ~= 0
        magnitude_ratio = fitted_peak / true_peak;
        magnitude_error = abs(fitted_peak - true_peak) / abs(true_peak);
    else
        magnitude_ratio = NaN;
        magnitude_error = NaN;
    end

    % store metrics
    metrics.(kernel_name) = struct();
    metrics.(kernel_name).correlation = corr_val;
    metrics.(kernel_name).peak_error_s = peak_error_s;
    metrics.(kernel_name).peak_error_bins = peak_error_bins;
    metrics.(kernel_name).fitted_peak = fitted_peak;
    metrics.(kernel_name).true_peak = true_peak;
    metrics.(kernel_name).magnitude_ratio = magnitude_ratio;
    metrics.(kernel_name).magnitude_error = magnitude_error;
end

% section print formatted table
fprintf('\n');
fprintf('========================================================================\n');
fprintf('                        KERNEL COMPARISON METRICS                       \n');
fprintf('========================================================================\n');
fprintf('\n');

fprintf('%-20s %12s %12s %12s %12s\n', ...
    'Kernel', 'Correlation', 'Peak Error', 'Magnitude', 'Mag Error');
fprintf('%-20s %12s %12s %12s %12s\n', ...
    '', '', '(bins)', 'Ratio', '(%)');
fprintf('------------------------------------------------------------------------\n');

for i = 1:length(kernel_names)
    kernel_name = kernel_names{i};

    if isfield(metrics, kernel_name)
        m = metrics.(kernel_name);

        % format output with pass/fail indicators
        corr_str = sprintf('%.4f', m.correlation);
        if m.correlation > 0.7
            corr_str = ['✓ ' corr_str];
        else
            corr_str = ['✗ ' corr_str];
        end

        peak_str = sprintf('%.2f', m.peak_error_bins);
        if m.peak_error_bins <= 2
            peak_str = ['✓ ' peak_str];
        else
            peak_str = ['✗ ' peak_str];
        end

        mag_str = sprintf('%.3f', m.magnitude_ratio);
        if m.magnitude_error < 0.4
            mag_str = ['✓ ' mag_str];
        else
            mag_str = ['✗ ' mag_str];
        end

        mag_err_str = sprintf('%.1f', m.magnitude_error * 100);

        fprintf('%-20s %12s %12s %12s %12s\n', ...
            kernel_name, corr_str, peak_str, mag_str, mag_err_str);
    else
        fprintf('%-20s %12s %12s %12s %12s\n', ...
            kernel_name, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end

fprintf('------------------------------------------------------------------------\n');
fprintf('\n');
fprintf('Criteria: Correlation > 0.7, Peak Error <= 2 bins, Magnitude Error < 40%%\n');
fprintf('========================================================================\n');
fprintf('\n');

end
