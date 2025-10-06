function tests = test_plot_outputs
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testPlotKernelsWritesFile(testCase)
% section kernels plot
% confirm kernel plotting saves an output image.
[kernels, outdir, cleanupObj] = makeKernelFixture(); %#ok<NASGU>
plot_kernels(kernels, outdir);
verifyFileExists(testCase, fullfile(outdir, 'kernels.pdf'));
end

function testPlotRateVsSpikesWritesFile(testCase)
% section rate plot
% confirm rate-vs-spikes plot saves an output image.
[~, outdir, cleanupObj] = makeKernelFixture(); %#ok<NASGU>
stim = struct('t', (0:0.1:0.5)', 'dt', 0.1);
y = [0; 1; 0; 1; 0; 1];
mu = [0.2; 0.8; 0.3; 0.9; 0.4; 0.7];
plot_rate_vs_spikes(stim, y, mu, outdir);
verifyFileExists(testCase, fullfile(outdir, 'rate_vs_spikes.pdf'));
end

function testPlotCvCurveWritesFile(testCase)
% section cv plot
% confirm cv curve plot saves an output image.
[~, outdir, cleanupObj] = makeKernelFixture(); %#ok<NASGU>
cvinfo = struct('lambdas', [0.1, 1, 10], ...
    'mean_nll', [1.0, 0.8, 1.1], ...
    'fold_nll', [1.1, 0.9, 1.2; 0.95, 0.75, 1.05; 1.05, 0.85, 1.15]);
plot_cv_curve(cvinfo, outdir);
verifyFileExists(testCase, fullfile(outdir, 'cv_curve.pdf'));
end

function testPlotPsthsWritesFile(testCase)
% section psth plot
% confirm the psth figure is written and handles produced subclasses.
[~, outdir, cleanupObj] = makeKernelFixture(); %#ok<NASGU>

sp = struct('spike_times', linspace(0, 2, 200));
heardEvents = struct('kind', 'perceived', 't_on', 0.3, 't_off', 0.35, 'label', "heard");
producedEvents = struct( ...
    'kind', {'produced', 'produced', 'produced'}, ...
    't_on', {0.05, 0.5, 0.7}, ...
    't_off', {0.10, 0.55, 0.75}, ...
    'label', {"p1", "p2", "p3"} ...
);

cfg = struct('heard_window_s', [-0.2, 0.4], 'produced_window_s', [-0.3, 0.5], 'dt', 0.01);

plot_psths(sp, heardEvents, producedEvents, cfg, outdir);
verifyFileExists(testCase, fullfile(outdir, 'psth.pdf'));
end

function verifyFileExists(testCase, filepath)
info = dir(filepath);
testCase.verifyTrue(~isempty(info) && info.bytes > 0);
end

function [kernels, outdir, cleanupObj] = makeKernelFixture()
outdir = tempname;
mkdir(outdir);
cleanupObj = onCleanup(@() rmdir(outdir, 's'));

kernels = struct();
kernels.intercept = 0.1;
kernels.heard_fields = {'heard_addressed', 'heard_overheard'};
kernels.heard_addressed = struct('weights', [0.1; 0.2], 'lag_times_s', [0; 0.01], 'mode', 'causal');
kernels.heard_overheard = struct('weights', [0.05; 0.01], 'lag_times_s', [0; 0.01], 'mode', 'causal');
kernels.produced_spontaneous = struct('weights', [-0.05; 0.05], 'lag_times_s', [-0.01; 0], 'mode', 'symmetric');
kernels.produced_after_heard = struct('weights', [-0.02; 0.03], 'lag_times_s', [-0.01; 0], 'mode', 'symmetric');
kernels.produced_after_produced = struct('weights', [0.04; -0.04], 'lag_times_s', [-0.01; 0], 'mode', 'symmetric');
kernels.spike_history = struct('weights', [0.3; 0.2; 0.1], 'lag_times_s', [0.01; 0.02; 0.03], 'mode', 'history');
kernels.states = struct('weights', [0.5; -0.2], 'names', {{'convo', 'spon'}}, 'convo', 0.5, 'spon', -0.2);
end
