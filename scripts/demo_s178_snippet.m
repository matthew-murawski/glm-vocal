function demo_s178_snippet()
% section overview
% run the single-neuron glm pipeline on the packaged 60 s s178 demo data.

[rootDir, ~] = resolve_paths();
dataDir = fullfile(rootDir, 'demos', 'data');

cfgPath = fullfile(rootDir, 'config', 'defaults.json');
spikePath = fullfile(dataDir, 'S178_demo_spike_times.mat');
heardPath = fullfile(dataDir, 'S178_test_heard.txt');
producedPath = fullfile(dataDir, 'S178_test_produced.txt');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(rootDir, 'results', ['demo_s178_' timestamp]);

fprintf('Running S178 demo → output: %s\n', outdir);
run_fit_single_neuron(cfgPath, spikePath, heardPath, producedPath, outdir);

fprintf('Finished. Inspect %s for artifacts.\n', outdir);
end

function [rootDir, srcDir] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcDir = genpath(fullfile(rootDir, 'src'));
addpath(srcDir);
end
