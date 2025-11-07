function results = run_fit_lfp_multichannel(cfgPath, lfpPath, heardPath, producedPath, outdir)
%RUN_FIT_LFP_MULTICHANNEL Fit GLM to multiple LFP channels.
%   results = RUN_FIT_LFP_MULTICHANNEL(cfgPath, lfpPath, heardPath, producedPath, outdir)
%   fits a Gaussian GLM to each LFP channel independently, using the same
%   regressors (heard/produced stimuli, conversational states) as the spike model.
%
%   Inputs:
%     cfgPath      - path to config JSON file (or struct)
%     lfpPath      - path to LFP MAT file
%     heardPath    - path to heard event labels (optional)
%     producedPath - path to produced event labels (optional)
%     outdir       - output directory for results (optional)
%
%   Returns:
%     results - array of structs (one per channel) with fitted models

% section setup
% coordinate the full multi-channel LFP GLM fitting pipeline and persist results.
if nargin < 1
    cfgPath = [];
end
if nargin < 2 || isempty(lfpPath)
    error('run_fit_lfp_multichannel:MissingLFPPath', 'LFP MAT file path is required.');
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

[rootDir, srcPath] = resolve_paths();
addpath(srcPath);

% load configuration
if isstruct(cfgPath)
    cfg = cfgPath;
else
    cfgFile = resolve_file(cfgPath, fullfile(rootDir, 'config', 'defaults.json'));
    cfg = jsondecode(fileread(cfgFile));
end

% set response type to LFP
cfg.response_type = 'lfp';

% resolve file paths
lfpFile = resolve_file(lfpPath);
heardFile = resolve_file(heardPath, '', true);
producedFile = resolve_file(producedPath, '', true);

% create output directory
if isempty(outdir)
    outdir = fullfile(rootDir, 'results_lfp');
end
ensure_dir(outdir);
plotDir = fullfile(outdir, 'plots');
ensure_dir(plotDir);

% configure predictor exclusions
if ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    cfg.exclude_predictors = {};
end

cfg.produced_split_mode = sanitize_produced_split_mode(cfg);
cfg.twitter_bout_window_s = sanitize_twitter_window(cfg);

% section load data
fprintf('Loading LFP data from %s...\n', lfpFile);
lfp_data = load_lfp(lfpFile);
fprintf('Loaded %d channels, %d samples, fs=%.1f Hz\n', ...
        lfp_data.n_channels, lfp_data.n_samples, lfp_data.fs);

% load event labels
heardEvents = repmat(empty_event_record(), 0, 1);
producedEvents = heardEvents;
if ~isempty(heardFile)
    heardEvents = load_labels(heardFile, 'perceived');
end
if ~isempty(producedFile)
    producedEvents = load_labels(producedFile, 'produced');
end
events = [heardEvents; producedEvents];

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
% create a dummy spike struct for build_timebase
dummy_sp = struct('spike_times', [], 'neuron_id', 'LFP', 'session_id', 'LFP');
stim = build_timebase(events, dummy_sp, cfg.dt);

% bin LFP data to timebase
fprintf('Binning LFP data to timebase (dt=%.3f s)...\n', cfg.dt);
lfp_binned = bin_lfp(lfp_data, stim, cfg);

% build streams and states (same as spike pipeline)
streams = build_streams(events, stim, cfg);
states = compute_states(events, stim, cfg.state);

% section fit each channel independently
n_channels = lfp_data.n_channels;
results = repmat(struct(), n_channels, 1);

fprintf('Fitting GLM to %d channels...\n', n_channels);
for ch = 1:n_channels
    fprintf('  Channel %d/%d (ID: %s)... ', ch, n_channels, num2str(lfp_data.channel_ids(ch)));

    lfp_ch = lfp_binned(:, ch);

    % assemble design matrix with LFP history
    Xd = assemble_design_matrix(streams, states, lfp_ch, cfg, stim);

    % build smoothness penalty matrix
    [D, Dmap] = smoothness_penalty(Xd.colmap, cfg);

    % get link function from config
    if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'link')
        link = cfg.lfp.link;
    else
        link = 'identity';
    end

    % cross-validate lambda
    fitfun = @(XdTrain, Dtrain, lambda) fit_glm_gauss(XdTrain, Dtrain, lambda, cfg.optimizer, link);
    [best_lambda, cvinfo] = crossval_blocked(fitfun, Xd, D, cfg.cv);

    % fit final model with best lambda
    [wmap, fitinfo] = fit_glm_gauss(Xd, D, best_lambda, cfg.optimizer, link);

    % predict LFP
    lfp_pred = predict_lfp(Xd.X, wmap.w, link);

    % compute metrics
    metricsOut = metrics_lfp(Xd.y, lfp_pred);
    fprintf('R²=%.3f, Corr=%.3f\n', metricsOut.r2, metricsOut.correlation);

    % unpack kernels
    kernels = unpack_params(wmap, Xd.colmap, cfg, stim);

    % store results
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
end

% section generate summary plots
fprintf('Generating summary plots...\n');

% plot R² across channels
fig_r2 = plot_lfp_r2_by_channel(results);
saveas(fig_r2, fullfile(plotDir, 'lfp_r2_by_channel.pdf'));
close(fig_r2);

% plot traces for a few representative channels
channels_to_plot = min(4, n_channels);
for ii = 1:channels_to_plot
    ch = ii;
    fig_trace = plot_lfp_traces(stim.t, results(ch).lfp_actual, results(ch).lfp_pred, ...
                                 results(ch).channel_id, results(ch).metrics);
    saveas(fig_trace, fullfile(plotDir, sprintf('lfp_traces_ch%d.pdf', ch)));
    close(fig_trace);
end

% section save results
fprintf('Saving results to %s...\n', outdir);
artifactPath = fullfile(outdir, 'fit_results_lfp.mat');
save(artifactPath, 'cfg', 'lfp_data', 'events', 'stim', 'lfp_binned', 'streams', 'states', 'results', '-v7.3');

% return consolidated results
result = struct();
result.cfg = cfg;
result.stim = stim;
result.results = results;

fprintf('Done! Results saved to: %s\n', artifactPath);

end
