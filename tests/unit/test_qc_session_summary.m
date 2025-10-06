function tests = test_qc_session_summary
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testSummaryArtifacts(testCase)
% section smoke test
% ensure the qc summary writes figure, text, and json outputs without error.
outdir = tempname;
mkdir(outdir);
cleanupObj = onCleanup(@() rmdir(outdir, 's'));

Xd = struct('X', sparse(rand(20, 5)), 'y', (0:19)' > 10);
wmap = struct('w', randn(5, 1));
rate = struct('stim', struct('t', (0:19)'), 'y', double((0:19)' > 5), 'mu', rand(20, 1));
rate.event_counts = struct('heard', 3, ...
    'heard_fields', {{'heard_addressed', 'heard_overheard'}}, ...
    'heard_by_field', struct('heard_addressed', 2, 'heard_overheard', 1), ...
    'produced_fields', {{'produced_spontaneous'}}, ...
    'produced', struct('produced_spontaneous', 2), ...
    'produced_any', 7);
cvinfo = struct('lambdas', [0.1, 1, 10], 'mean_nll', [1.0, 0.8, 0.95]);

summary = qc_session_summary(Xd, wmap, rate, cvinfo, struct(), struct(), outdir);

verifyFile(testCase, summary.paths.figure);
verifyFile(testCase, summary.paths.text);
verifyFile(testCase, summary.paths.json);

testCase.verifyEqual(summary.stats.event_counts.heard, rate.event_counts.heard);
testCase.verifyEqual(summary.stats.event_counts.produced.produced_spontaneous, 2);

clear cleanupObj; %#ok<CLCLR>
end

function testPermutationStatsIncluded(testCase)
% section permutation summary
% validate that permutation test results are copied into the qc summary structure.
outdir = tempname;
mkdir(outdir);
cleanupObj = onCleanup(@() rmdir(outdir, 's'));

Xd = struct('X', sparse(rand(10, 3)));
wmap = struct('w', randn(3, 1));
rate = struct();
cvinfo = struct();
kernels = struct();
ptest = struct('heard_addressed', struct('p_value', 0.02, 'ci_lower', [-0.1, -0.05], 'ci_upper', [0.05, 0.1]));

summary = qc_session_summary(Xd, wmap, rate, cvinfo, kernels, ptest, outdir);

permStats = summary.stats.permutation;
testCase.verifyEqual(permStats.n_kernels, 1);
testCase.verifyEqual(permStats.n_significant, 1);
testCase.verifyTrue(isfield(permStats.kernels, 'heard_addressed'));
testCase.verifyLessThan(permStats.kernels.heard_addressed.p_value, permStats.threshold);

clear cleanupObj; %#ok<CLCLR>
end

function verifyFile(testCase, filepath)
info = dir(filepath);
testCase.verifyTrue(~isempty(info) && info.bytes > 0);
end
