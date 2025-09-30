function demo_s177()
% section overview
% run the single-neuron glm pipeline on all of s177.

[rootDir, ~] = resolve_paths();
dataDir = '/Users/matt/Documents/GitHub/vocalization/data/Label Files/S177';

cfgPath = fullfile(rootDir, 'config', 'defaults.json');
spikePath = fullfile(dataDir, 'M93A_S177_spike_times.mat');
heardPath = fullfile(dataDir, 'M93A_S177_heard.txt');
producedPath = fullfile(dataDir, 'M93A_S177_produced.txt');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(rootDir, 'results', ['S177_' timestamp]);

fprintf('Running S177 → output: %s\n', outdir);
run_fit_single_neuron(cfgPath, spikePath, heardPath, producedPath, outdir);

fprintf('Finished. Inspect %s for artifacts.\n', outdir);
end

function [rootDir, srcDir] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcDir = genpath(fullfile(rootDir, 'src'));
addpath(srcDir);
end