function plot_rate_vs_spikes(stim, y, mu, outdir)
% section setup
% produce smooth traces (no vertical lines) by gaussian-smoothing both series and resampling at 100 ms.
if nargin < 4 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% basic checks
t = stim.t(:);
if numel(t) ~= numel(y) || numel(y) ~= numel(mu)
    error('plot_rate_vs_spikes:SizeMismatch', 'stim, y, and mu must align.');
end
if ~isfield(stim, 'dt') || ~isscalar(stim.dt) || ~isfinite(stim.dt) || stim.dt <= 0
    error('plot_rate_vs_spikes:BadDt', 'stim.dt must be a positive scalar.');
end

y  = y(:);                 % counts per native bin
mu = mu(:);                % predicted (units specified below)

% section units policy
% set to false if mu is expected spikes per native bin instead of spikes/s.
mu_is_rate = false;         % assumption: mu is spikes/s

% section smoothing configuration
% gaussian kernel: sigma = 50 ms, length = ~6*sigma (odd), native grid at stim.dt
dt              = stim.dt;
sigma_sec       = 0.050;                                % 50 ms
sigma_samp      = max(1, round(sigma_sec / dt));
halfwidth_samp  = 3 * sigma_samp;                       % ~3 sigma on each side
x               = -halfwidth_samp:halfwidth_samp;
g               = exp(-0.5 * (x / sigma_samp).^2);
g               = g / sum(g);                           % normalize to 1

% section smooth to rate
% observed: counts/bin -> smooth -> divide by dt => spikes/s
y_rate_native = conv(y, g, 'same') ./ dt;               % spikes/s, smooth

% predicted: match bandwidth and units
if mu_is_rate
    mu_rate_native = conv(mu, g, 'same');               % already spikes/s
else
    mu_rate_native = conv(mu, g, 'same') ./ dt;         % per-bin -> spikes/s
end

% section resample to 100 ms for a clean trace
% choose centers aligned to kernel midpoint to avoid phase lag
target_step_sec = 0.1;                                  % 100 ms
k               = max(1, round(target_step_sec / dt));
start_idx       = ceil(numel(g)/2);                     % center-aligned sampling
idx             = start_idx:k:(numel(t) - (ceil(numel(g)/2)-1));
t_ds            = t(idx);
y_rate_ds       = y_rate_native(idx);
mu_rate_ds      = mu_rate_native(idx);

% section plot
fig = figure('Visible', 'off');
try
    % use line plots only; no stems/bars; moderate linewidth
    plot(t_ds, y_rate_ds, 'LineWidth', 1.4, 'DisplayName', 'observed (spikes/s, 100 ms σ=50 ms)');
    hold on
    plot(t_ds, mu_rate_ds, 'LineWidth', 1.4, 'DisplayName', 'predicted (spikes/s, 100 ms σ=50 ms)');
    hold off
    xlabel('time (s)');
    ylabel('spikes/s');
    title('rate vs spikes — gaussian smoothed and resampled (no vertical lines)');
    legend('Location', 'best');
    grid on;

    % export; png avoids hairline aliasing that can look like "pickets" in some pdf viewers
    exportgraphics(fig, fullfile(outdir, 'rate_vs_spikes.png'), 'Resolution', 300);
    exportgraphics(fig, fullfile(outdir, 'rate_vs_spikes.pdf'), 'ContentType', 'vector');
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end