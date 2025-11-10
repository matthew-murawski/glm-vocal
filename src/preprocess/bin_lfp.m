function lfp_binned = bin_lfp(lfp, stim, cfg)
%BIN_LFP Resample continuous LFP data onto the stimulus time grid.
%   lfp_binned = BIN_LFP(lfp, stim, cfg) resamples LFP data to match the
%   stimulus timebase, with optional filtering and normalization.
%
%   Inputs:
%     lfp  - struct from load_lfp with fields:
%            .data (n_samples × n_channels)
%            .fs (sampling rate in Hz)
%            .t0 (start time in seconds)
%     stim - struct with fields:
%            .t (time bin centers, column vector)
%            .dt (bin width in seconds)
%     cfg  - config struct with optional fields:
%            .lfp.bandpass_hz (2-element vector [low high] or empty)
%            .lfp.normalize ('zscore', 'minmax', or 'none')
%
%   Returns:
%     lfp_binned - (n_bins × n_channels) matrix aligned to stim.t

% validate inputs
if nargin < 3
    cfg = struct();
end

% construct the original time vector for the LFP data
t_lfp = lfp.t0 + (0:(lfp.n_samples-1))' / lfp.fs;

% extract configuration options
if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'bandpass_hz') && ~isempty(cfg.lfp.bandpass_hz)
    bandpass_hz = cfg.lfp.bandpass_hz;
else
    bandpass_hz = [];
end

if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'normalize')
    normalize_mode = cfg.lfp.normalize;
else
    normalize_mode = 'zscore';  % default
end

% optional bandpass filtering
if ~isempty(bandpass_hz)
    if numel(bandpass_hz) ~= 2 || bandpass_hz(1) >= bandpass_hz(2)
        error('glm:InvalidBandpass', 'Bandpass must be a 2-element vector [low high] with low < high.');
    end

    % design butterworth bandpass filter
    filter_order = 4;
    nyquist = lfp.fs / 2;
    Wn = bandpass_hz / nyquist;

    % handle edge case where high cutoff exceeds Nyquist
    if Wn(2) >= 1
        warning('glm:BandpassTooHigh', 'High cutoff (%g Hz) >= Nyquist (%g Hz). Using highpass only.', ...
                bandpass_hz(2), nyquist);
        [b, a] = butter(filter_order, Wn(1), 'high');
    else
        [b, a] = butter(filter_order, Wn, 'bandpass');
    end

    % apply filter to each channel
    lfp_filtered = zeros(size(lfp.data));
    for ch = 1:lfp.n_channels
        lfp_filtered(:, ch) = filtfilt(b, a, lfp.data(:, ch));
    end
else
    lfp_filtered = lfp.data;
end

% extract band power if requested
if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'extract_band_power') && cfg.lfp.extract_band_power
    fprintf('  Extracting band power via Hilbert transform...\n');

    % get smoothing window if specified
    if isfield(cfg.lfp, 'power_smoothing_ms') && ~isempty(cfg.lfp.power_smoothing_ms)
        smoothing_ms = cfg.lfp.power_smoothing_ms;
    else
        smoothing_ms = 0;  % no smoothing
    end

    % extract power for each channel
    lfp_power = zeros(size(lfp_filtered));
    for ch = 1:lfp.n_channels
        % apply Hilbert transform to get analytic signal
        analytic_signal = hilbert(lfp_filtered(:, ch));

        % extract amplitude envelope (instantaneous amplitude)
        envelope = abs(analytic_signal);

        % optional smoothing of the envelope
        if smoothing_ms > 0
            % convert smoothing window from ms to samples
            smoothing_samples = round(smoothing_ms * lfp.fs / 1000);

            % use moving average filter for smoothing
            if smoothing_samples > 1
                smooth_kernel = ones(smoothing_samples, 1) / smoothing_samples;
                envelope = filtfilt(smooth_kernel, 1, envelope);
            end
        end

        lfp_power(:, ch) = envelope;
    end

    % use power envelope for subsequent processing
    lfp_filtered = lfp_power;
end

% resample to stimulus grid using interpolation
n_bins = numel(stim.t);
lfp_binned = zeros(n_bins, lfp.n_channels);

for ch = 1:lfp.n_channels
    % use interp1 with linear interpolation (fast and robust)
    % extrapolation: use 0 for times outside the LFP recording window
    lfp_binned(:, ch) = interp1(t_lfp, lfp_filtered(:, ch), stim.t, 'linear', 0);
end

% optional normalization per channel
switch normalize_mode
    case 'zscore'
        % z-score: subtract mean, divide by std
        for ch = 1:lfp.n_channels
            mu = mean(lfp_binned(:, ch));
            sigma = std(lfp_binned(:, ch));
            if sigma > 0
                lfp_binned(:, ch) = (lfp_binned(:, ch) - mu) / sigma;
            else
                warning('glm:ZeroStd', 'Channel %d has zero std; skipping normalization.', ch);
            end
        end

    case 'minmax'
        % min-max normalization to [0, 1]
        for ch = 1:lfp.n_channels
            minval = min(lfp_binned(:, ch));
            maxval = max(lfp_binned(:, ch));
            if maxval > minval
                lfp_binned(:, ch) = (lfp_binned(:, ch) - minval) / (maxval - minval);
            else
                warning('glm:ConstantSignal', 'Channel %d is constant; skipping normalization.', ch);
            end
        end

    case 'none'
        % no normalization

    otherwise
        error('glm:InvalidNormalize', 'Unknown normalization mode: %s', normalize_mode);
end

end
