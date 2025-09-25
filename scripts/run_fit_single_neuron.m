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

cfgFile = resolve_file(cfgPath, fullfile(rootDir, 'config', 'defaults.json'));
spikeFile = resolve_file(spikePath);
heardFile = resolve_file(heardPath, '', true);
producedFile = resolve_file(producedPath, '', true);

if isempty(outdir)
    outdir = fullfile(rootDir, 'results');
end
ensure_dir(outdir);
plotDir = fullfile(outdir, 'plots');
ensure_dir(plotDir);

cfg = jsondecode(fileread(cfgFile));
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

stim = build_timebase(events, sp, cfg.dt);
sps = bin_spikes(sp.spike_times, stim);
streams = build_streams(events, stim);
states = compute_states(events, stim, cfg.state);

Xd = assemble_design_matrix(streams, states, sps, cfg, stim);
[D, Dmap] = smoothness_penalty(Xd.colmap, cfg);

fitfun = @(XdTrain, Dtrain, lambda) fit_glm_map(XdTrain, Dtrain, lambda, cfg.optimizer);
[best_lambda, cvinfo] = crossval_blocked(@(XdTrain, Dtrain, lambda) fitfun(XdTrain, Dtrain, lambda), Xd, D, cfg.cv);
[wmap, fitinfo] = fit_glm_map(Xd, D, best_lambda, cfg.optimizer);

mu = predict_rate(Xd.X, wmap.w);
kernels = unpack_params(wmap, Xd.colmap, cfg, stim);
metricsOut = metrics(Xd.y, mu);
rate = struct('stim', stim, 'y', Xd.y, 'mu', mu, 'metrics', metricsOut);

plot_design_matrix(Xd.X, Xd.colmap, struct('output_path', fullfile(plotDir, 'design_matrix.png')));
plot_kernels(kernels, plotDir);
plot_rate_vs_spikes(stim, Xd.y, mu, plotDir);
plot_cv_curve(cvinfo, plotDir);

summary = qc_session_summary(Xd, wmap, rate, cvinfo, outdir);

artifactPath = fullfile(outdir, 'fit_results.mat');
save(artifactPath, 'cfg', 'sp', 'events', 'stim', 'sps', 'streams', 'states', 'Xd', 'D', 'Dmap', 'best_lambda', 'cvinfo', 'wmap', 'fitinfo', 'kernels', 'metricsOut', 'summary');

if nargout > 0
    result = struct('cfg', cfg, 'sp', sp, 'events', events, 'Xd', Xd, 'wmap', wmap, ...
        'cvinfo', cvinfo, 'fitinfo', fitinfo, 'kernels', kernels, 'metrics', metricsOut, ...
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

function rec = empty_event_record()
rec = struct('kind', '', 't_on', [], 't_off', [], 'label', "");
end
