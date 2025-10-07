function demo_s177()
% section overview
% run the single-neuron glm pipeline on all of s177.

[rootDir, ~] = resolve_paths();
dataDir = '/Users/matt/Documents/GitHub/vocalization/data/Label Files/S177';
cfgPath = fullfile(rootDir, 'config', 'defaults.json');
spikePath = fullfile(dataDir, 'M93A_S177_spike_times_ch32.mat');
heardPath = fullfile(dataDir, 'M93A_S177_heard.txt');
producedPath = fullfile(dataDir, 'M93A_S177_produced.txt');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(rootDir, 'results', ['S177_' timestamp]);

cfg = jsondecode(fileread(cfgPath));
if ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    cfg.exclude_predictors = {'states', 'spike_history'};
end
cfg.produced_split_mode = 'call_type';

fprintf('Running S177 → output: %s\n', outdir);
result = run_fit_single_neuron(cfg, spikePath, heardPath, producedPath, outdir);
fprintf('Finished. Inspect %s for artifacts.\n', outdir);

% section produced-call preview
plot_event_prediction_preview(result); % Uses all defaults

% plot_event_prediction_preview(result, 'EventType', 'perceived', 'EventIndex', 5);

end

function [rootDir, srcDir] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcDir = genpath(fullfile(rootDir, 'src'));
addpath(srcDir);
end