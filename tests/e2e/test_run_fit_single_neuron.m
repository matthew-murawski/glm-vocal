function tests = test_run_fit_single_neuron
% section registration
% integration smoke test for the single-neuron fitting pipeline.
tests = functiontests(localfunctions);
end

function testSyntheticPipelineProducesArtifacts(testCase)
% section synthetic pipeline
% run the full pipeline on synthetic data and confirm key artifacts exist.
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'scripts'));

workDir = tempname;
mkdir(workDir);
cleanupObj = onCleanup(@() rmdir(workDir, 's'));

cfgStruct = struct();
cfgStruct.dt = 0.05;
cfgStruct.heard_window_s = [0, 0.2];
cfgStruct.produced_window_s = [-0.2, 0.2];
cfgStruct.history_window_s = [0.05, 0.3];
cfgStruct.cv = struct('k', 2, 'lambdas', [0.5]);
cfgStruct.state = struct('response_window_s', 2.0, 'max_seq_gap_s', 2.0);
cfgStruct.optimizer = struct('tol_fun', 1e-6, 'tol_grad', 1e-5, 'max_iter', 200, 'use_fminunc', false);
cfgStruct.seed = 123;
cfgPath = fullfile(workDir, 'config.json');
fidCfg = fopen(cfgPath, 'w');
fwrite(fidCfg, jsonencode(cfgStruct), 'char');
fclose(fidCfg);

% synthetic spikes (seconds)
spike_times = sort(rand(100, 1) * 2);
neuron_id = "demo";
session_id = "synthetic";
spikeFile = fullfile(workDir, 'spikes.mat');
save(spikeFile, 'spike_times', 'neuron_id', 'session_id');

% simple heard/produced label txt files (start stop label)
heardFile = fullfile(workDir, 'heard.txt');
writeLabelFile(heardFile, [0.2, 0.4; 0.9, 1.0], 'perceived');
producedFile = fullfile(workDir, 'produced.txt');
writeLabelFile(producedFile, [0.5, 0.6; 1.2, 1.3], 'produced');

outdir = fullfile(workDir, 'results');
result = run_fit_single_neuron(cfgPath, spikeFile, heardFile, producedFile, outdir);

% verify artifacts
verifyFile(testCase, fullfile(outdir, 'fit_results.mat'));
verifyFile(testCase, fullfile(outdir, 'qc_summary.png'));
verifyFile(testCase, fullfile(outdir, 'qc_summary.txt'));
verifyFile(testCase, fullfile(outdir, 'qc_summary.json'));
verifyFile(testCase, fullfile(outdir, 'plots', 'kernels.pdf'));
verifyFile(testCase, fullfile(outdir, 'plots', 'rate_vs_spikes.pdf'));
verifyFile(testCase, fullfile(outdir, 'plots', 'cv_curve.pdf'));

% basic sanity checks on outputs
fitData = load(fullfile(outdir, 'fit_results.mat'));
testCase.verifyTrue(isfield(fitData, 'wmap'));
testCase.verifyTrue(isfield(result, 'metrics'));
testCase.verifyGreaterThan(result.metrics.pseudoR2, -Inf);
end

function writeLabelFile(pathStr, intervals, label)
fid = fopen(pathStr, 'w');
for ii = 1:size(intervals, 1)
    fprintf(fid, '%.6f\t%.6f\t%s\n', intervals(ii, 1), intervals(ii, 2), label);
end
fclose(fid);
end

function verifyFile(testCase, filepath)
info = dir(filepath);
testCase.verifyTrue(~isempty(info) && info.bytes > 0);
end
