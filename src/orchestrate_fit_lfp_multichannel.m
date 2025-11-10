function results = orchestrate_fit_lfp_multichannel(cfgPath, lfpPath, heardPath, producedPath, outdir, P)

% section setup
if nargin < 1
    cfgPath = [];
end
if nargin < 2 || isempty(lfpPath)
    error('orchestrate_fit_lfp_multichannel:MissingLFPPath', 'LFP MAT file path is required.');
end
if nargin < 3
    heardPath = [];
end
if nargin < 4
    producedPath = [];
end
if nargin < 5
    outdir = [];
end

rootDir = fullfile(P.github_path, 'glm-vocal');
srcPath = fullfile(P.github_path, 'glm-vocal/src');

addpath(srcPath);

% load configuration
if isstruct(cfgPath)
    cfg = cfgPath;
elseif ~isempty(cfgPath)
    cfg = jsondecode(fileread(cfgPath));
else
    cfgFile = fullfile(rootDir, 'config', 'defaults.json');
    cfg = jsondecode(fileread(cfgFile));
end

% set response type to LFP
cfg.response_type = 'lfp';

% load event labels
heardEvents = repmat(empty_event_record(), 0, 1);
producedEvents = heardEvents;
if ~isempty(heardPath) && isfile(heardPath)
    heardEvents = load_labels(heardPath, 'perceived');
end
if ~isempty(producedPath) && isfile(producedPath)
    producedEvents = load_labels(producedPath, 'produced');
end
events = [heardEvents; producedEvents];

% create output directory
if isempty(outdir)
    outdir = fullfile(rootDir, 'results_lfp');
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
plotDir = fullfile(outdir, 'plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% configure predictor exclusions
if ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    cfg.exclude_predictors = {};
end

cfg.produced_split_mode = sanitize_produced_split_mode(cfg);
cfg.twitter_bout_window_s = sanitize_twitter_window(cfg);

% section load data
fprintf('Loading LFP data from %s...\n', lfpPath);
lfp_data = load_lfp(lfpPath);
fprintf('Loaded %d channels, %d samples, fs=%.1f Hz\n', ...
        lfp_data.n_channels, lfp_data.n_samples, lfp_data.fs);

% consolidate twitter bouts
events = consolidate_twitter_bouts(events, cfg.twitter_bout_window_s);
if ~isempty(events)
    kindCells = cellfun(@char, {events.kind}, 'UniformOutput', false);
    heardEvents = events(strcmp(kindCells, 'perceived'));
    producedEvents = events(strcmp(kindCells, 'produced'));
else
    heardEvents = repmat(empty_event_record(), 0, 1);
    producedEvents = heardEvents;
end

% section build shared preprocessing structures
fprintf('Building timebase and streams...\n');
dummy_sp = struct('spike_times', [], 'neuron_id', 'LFP', 'session_id', 'LFP');
stim = build_timebase(events, dummy_sp, cfg.dt);

fprintf('Binning LFP data to timebase (dt=%.3f s)...\n', cfg.dt);
lfp_binned = bin_lfp(lfp_data, stim, cfg);

streams = build_streams(events, stim, cfg);
states = compute_states(events, stim, cfg.state);

% section fit each channel independently
n_channels = lfp_data.n_channels;
results = repmat(struct(), n_channels, 1);

fprintf('Fitting GLM to %d channels...\n', n_channels);
for ch = 1:n_channels
    fprintf('  Channel %d/%d (ID: %s)... ', ch, n_channels, num2str(lfp_data.channel_ids(ch)));

    lfp_ch = lfp_binned(:, ch);

    Xd = assemble_design_matrix(streams, states, lfp_ch, cfg, stim);

    [D, Dmap] = smoothness_penalty(Xd.colmap, cfg);

    if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'link')
        link = cfg.lfp.link;
    else
        link = 'identity';
    end

    fitfun = @(XdTrain, Dtrain, lambda) fit_glm_gauss(XdTrain, Dtrain, lambda, cfg.optimizer, link);
    [best_lambda, cvinfo] = crossval_blocked(fitfun, Xd, D, cfg.cv);

    [wmap, fitinfo] = fit_glm_gauss(Xd, D, best_lambda, cfg.optimizer, link);

    lfp_pred = predict_lfp(Xd.X, wmap.w, link);

    metricsOut = metrics_lfp(Xd.y, lfp_pred);
    fprintf('R²=%.3f, Corr=%.3f\n', metricsOut.r2, metricsOut.correlation);

    kernels = unpack_params(wmap, Xd.colmap, cfg, stim);

    % perform permutation test for statistical significance
    ptest = perform_permutation_test_lfp(kernels, Xd, D, best_lambda, cfg, stim);

    results(ch).channel_id = lfp_data.channel_ids(ch);
    results(ch).Xd = Xd;
    results(ch).D = D;
    results(ch).Dmap = Dmap;
    results(ch).best_lambda = best_lambda;
    results(ch).cvinfo = cvinfo;
    results(ch).wmap = wmap;
    results(ch).fitinfo = fitinfo;
    results(ch).lfp_pred = lfp_pred;
    results(ch).lfp_actual = Xd.y;
    results(ch).metrics = metricsOut;
    results(ch).kernels = kernels;
    results(ch).ptest = ptest;
end

% section generate summary plots
fprintf('Generating summary plots...\n');

fig_r2 = plot_lfp_r2_by_channel(results);
saveas(fig_r2, fullfile(plotDir, 'lfp_r2_by_channel.pdf'));
close(fig_r2);

channels_to_plot = min(4, n_channels);
for ii = 1:channels_to_plot
    ch = ii;
    fig_trace = plot_lfp_traces(stim.t, results(ch).lfp_actual, results(ch).lfp_pred, ...
                                 results(ch).channel_id, results(ch).metrics);
    saveas(fig_trace, fullfile(plotDir, sprintf('lfp_traces_ch%d.pdf', ch)));
    close(fig_trace);
end

% plot fitted kernels for each channel
fprintf('Generating kernel plots...\n');
for ii = 1:n_channels
    ch = ii;

    % create detailed kernel plots in subdirectory
    kernel_outdir = fullfile(plotDir, sprintf('channel_%d', ch));
    if ~exist(kernel_outdir, 'dir')
        mkdir(kernel_outdir);
    end
    plot_kernels(results(ch).kernels, results(ch).ptest, kernel_outdir);

    % create summary kernel plot
    fig_kernels = plot_lfp_kernels_summary(results(ch).kernels, ...
                                           results(ch).channel_id, ...
                                           results(ch).metrics);
    saveas(fig_kernels, fullfile(plotDir, sprintf('lfp_kernels_summary_ch%d.pdf', ch)));
    close(fig_kernels);

    fprintf('  Kernels for channel %d saved\n', ch);
end

% section save results
fprintf('Saving results to %s...\n', outdir);
artifactPath = fullfile(outdir, 'fit_results_lfp.mat');
save(artifactPath, 'cfg', 'lfp_data', 'events', 'stim', 'lfp_binned', 'streams', 'states', 'results', '-v7.3');

result = struct();
result.cfg = cfg;
result.stim = stim;
result.results = results;

fprintf('Done! Results saved to: %s\n', artifactPath);

end

function mode = sanitize_produced_split_mode(cfg)
if ~isstruct(cfg) || ~isfield(cfg, 'produced_split_mode') || isempty(cfg.produced_split_mode)
    mode = 'context';
    return;
end
raw = string(cfg.produced_split_mode);
if numel(raw) ~= 1
    error('orchestrate_fit_lfp_multichannel:InvalidSplitMode', 'produced_split_mode must be a scalar string.');
end
mode = lower(strtrim(raw));
if ~(mode == "context" || mode == "call_type")
    error('orchestrate_fit_lfp_multichannel:InvalidSplitMode', 'produced_split_mode must be ''context'' or ''call_type''.');
end
mode = char(mode);
end

function window = sanitize_twitter_window(cfg)
defaultWindow = 1.5;
if ~isstruct(cfg) || ~isfield(cfg, 'twitter_bout_window_s') || isempty(cfg.twitter_bout_window_s)
    window = defaultWindow;
    return
end

raw = double(cfg.twitter_bout_window_s);
if ~isscalar(raw) || ~isfinite(raw) || raw <= 0
    error('orchestrate_fit_lfp_multichannel:InvalidTwitterWindow', 'twitter_bout_window_s must be a positive finite scalar.');
end
window = raw;
end