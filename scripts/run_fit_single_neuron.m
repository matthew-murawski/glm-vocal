function result = run_fit_single_neuron(cfgPath, spikePath, heardPath, producedPath, outdir)
% section setup
% coordinate the full single-neuron glm fitting pipeline and persist qc artifacts.
if nargin < 1
    cfgPath = [];
end
if nargin < 2 || isempty(spikePath)
    error('run_fit_single_neuron:MissingSpikePath', 'Spike MAT file path is required.');
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

if isstruct(cfgPath)
    cfg = cfgPath;
else
    cfgFile = resolve_file(cfgPath, fullfile(rootDir, 'config', 'defaults.json'));
    cfg = jsondecode(fileread(cfgFile));
end
spikeFile = resolve_file(spikePath);
heardFile = resolve_file(heardPath, '', true);
producedFile = resolve_file(producedPath, '', true);

if isempty(outdir)
    outdir = fullfile(rootDir, 'results');
end
ensure_dir(outdir);
plotDir = fullfile(outdir, 'plots');
ensure_dir(plotDir);

if ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    cfg.exclude_predictors = {};
end

cfg.produced_split_mode = sanitize_produced_split_mode(cfg);
cfg.twitter_bout_window_s = sanitize_twitter_window(cfg);
sp = load_spikes(spikeFile);

heardEvents = repmat(empty_event_record(), 0, 1);
producedEvents = heardEvents;
if ~isempty(heardFile)
    heardEvents = load_labels(heardFile, 'perceived');
end
if ~isempty(producedFile)
    producedEvents = load_labels(producedFile, 'produced');
end
events = [heardEvents; producedEvents];

events = consolidate_twitter_bouts(events, cfg.twitter_bout_window_s);
if ~isempty(events)
    kindCells = cellfun(@char, {events.kind}, 'UniformOutput', false);
    heardEvents = events(strcmp(kindCells, 'perceived'));
    producedEvents = events(strcmp(kindCells, 'produced'));
else
    heardEvents = repmat(empty_event_record(), 0, 1);
    producedEvents = heardEvents;
end

stim = build_timebase(events, sp, cfg.dt);
sps = bin_spikes(sp.spike_times, stim);
streams = build_streams(events, stim, cfg);
states = compute_states(events, stim, cfg.state);

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);
[D, Dmap] = smoothness_penalty(Xd.colmap, cfg);

fitfun = @(XdTrain, Dtrain, lambda) fit_glm_map(XdTrain, Dtrain, lambda, cfg.optimizer);
[best_lambda, cvinfo] = crossval_blocked(@(XdTrain, Dtrain, lambda) fitfun(XdTrain, Dtrain, lambda), Xd, D, cfg.cv);
[wmap, fitinfo] = fit_glm_map(Xd, D, best_lambda, cfg.optimizer);

mu = predict_rate(Xd.X, wmap.w);
kernels = unpack_params(wmap, Xd.colmap, cfg, stim);
ptest = perform_permutation_test(kernels, Xd, D, best_lambda, cfg, stim);
metricsOut = metrics(Xd.y, mu);
eventCounts = summarize_event_counts(streams);
rate = struct('stim', stim, 'y', Xd.y, 'mu', mu, 'metrics', metricsOut, 'event_counts', eventCounts);

plot_design_matrix(Xd.X, Xd.colmap, struct('output_path', fullfile(plotDir, 'design_matrix.pdf')));
plot_kernels(kernels, ptest, plotDir);
plot_rate_vs_spikes(stim, Xd.y, mu, plotDir);
plot_cv_curve(cvinfo, plotDir);
plot_psths(sp, heardEvents, producedEvents, cfg, plotDir);

summary = qc_session_summary(Xd, wmap, rate, cvinfo, kernels, ptest, outdir);

artifactPath = fullfile(outdir, 'fit_results.mat');
save(artifactPath, 'cfg', 'sp', 'events', 'stim', 'sps', 'streams', 'states', 'Xd', 'D', 'Dmap', 'best_lambda', 'cvinfo', 'wmap', 'fitinfo', 'kernels', 'ptest', 'metricsOut', 'summary');

if nargout > 0
    result = struct('cfg', cfg, 'sp', sp, 'events', events, 'Xd', Xd, 'wmap', wmap, ...
        'cvinfo', cvinfo, 'fitinfo', fitinfo, 'kernels', kernels, 'ptest', ptest, 'metrics', metricsOut, ...
        'summary', summary, 'paths', struct('outdir', outdir, 'artifact', artifactPath, 'plots', plotDir));
end
end

function [rootDir, srcPath] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcPath = genpath(fullfile(rootDir, 'src'));
end

function pathOut = resolve_file(pathIn, defaultPath, allowEmpty)
if nargin < 2
    defaultPath = '';
end
if nargin < 3
    allowEmpty = false;
end

if isempty(pathIn)
    pathOut = defaultPath;
else
    if isstring(pathIn)
        pathIn = char(pathIn);
    end
    if exist(pathIn, 'file') == 2
        pathOut = pathIn;
    else
        error('run_fit_single_neuron:MissingFile', 'File not found: %s', pathIn);
    end
end

if isempty(pathOut)
    if allowEmpty
        return
    else
        error('run_fit_single_neuron:MissingFile', 'Required file path not provided.');
    end
end
end

function ensure_dir(pathStr)
if exist(pathStr, 'dir') ~= 7
    mkdir(pathStr);
end
end

function counts = summarize_event_counts(streams)
% section event counts
% collect heard counts plus produced-call counts for each configured category.
counts = struct('heard', NaN, 'produced_fields', {{}}, 'produced', struct(), 'produced_any', NaN);

if ~isstruct(streams)
    return
end

if isfield(streams, 'heard_any')
    counts.heard = sum(double(streams.heard_any(:)) > 0);
end

if isfield(streams, 'produced_any')
    counts.produced_any = sum(double(streams.produced_any(:)) > 0);
end

producedFields = {};
if isfield(streams, 'produced_fields') && ~isempty(streams.produced_fields)
    producedFields = cellstr(streams.produced_fields(:));
else
    defaults = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
    producedFields = defaults(isfield(streams, defaults));
end

counts.produced_fields = producedFields;
for ii = 1:numel(producedFields)
    fieldName = producedFields{ii};
    if isfield(streams, fieldName)
        counts.produced.(fieldName) = sum(double(streams.(fieldName)(:)) > 0);
    else
        counts.produced.(fieldName) = NaN;
    end
end
end

function mode = sanitize_produced_split_mode(cfg)
if ~isstruct(cfg) || ~isfield(cfg, 'produced_split_mode') || isempty(cfg.produced_split_mode)
    mode = 'context';
    return
end

raw = string(cfg.produced_split_mode);
if numel(raw) ~= 1
    error('run_fit_single_neuron:InvalidSplitMode', 'produced_split_mode must be a scalar string.');
end

modeLower = lower(strtrim(raw));
validModes = ["context", "call_type"];
if ~any(modeLower == validModes)
    error('run_fit_single_neuron:InvalidSplitMode', 'produced_split_mode must be ''context'' or ''call_type''.');
end

mode = char(modeLower);
end

function window = sanitize_twitter_window(cfg)
defaultWindow = 1.5;
if ~isstruct(cfg) || ~isfield(cfg, 'twitter_bout_window_s') || isempty(cfg.twitter_bout_window_s)
    window = defaultWindow;
    return
end

raw = double(cfg.twitter_bout_window_s);
if ~isscalar(raw) || ~isfinite(raw) || raw <= 0
    error('run_fit_single_neuron:InvalidTwitterWindow', 'twitter_bout_window_s must be a positive finite scalar.');
end
window = raw;
end

function rec = empty_event_record()
rec = struct('kind', '', 't_on', [], 't_off', [], 'label', "");
end
