function demo_s178()
% section overview
% run the single-neuron glm pipeline on the packaged 60 s s178 demo data.

[rootDir, ~] = resolve_paths();
dataDir = '/Users/matt/Documents/GitHub/vocalization/data/demos/';

cfgPath = fullfile(rootDir, 'config', 'defaults.json');
spikePath = fullfile(dataDir, 'M93A_S178_spike_times.mat');
heardPath = fullfile(dataDir, 'M93A_S178_heard.txt');
producedPath = fullfile(dataDir, 'M93A_S178_produced.txt');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(rootDir, 'results', ['S178_' timestamp]);

fprintf('Running S178 → output: %s\n', outdir);
run_fit_single_neuron(cfgPath, spikePath, heardPath, producedPath, outdir);

fprintf('Finished. Inspect %s for artifacts.\n', outdir);
end

function [rootDir, srcDir] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcDir = genpath(fullfile(rootDir, 'src'));
addpath(srcDir);
end