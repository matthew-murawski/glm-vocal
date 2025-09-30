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
rate.event_counts = struct('heard', 3, 'produced_spontaneous', 2, ...
    'produced_after_heard', 4, 'produced_after_produced', 1);
cvinfo = struct('lambdas', [0.1, 1, 10], 'mean_nll', [1.0, 0.8, 0.95]);

summary = qc_session_summary(Xd, wmap, rate, cvinfo, struct(), outdir);

verifyFile(testCase, summary.paths.figure);
verifyFile(testCase, summary.paths.text);
verifyFile(testCase, summary.paths.json);

testCase.verifyEqual(summary.stats.event_counts.heard, rate.event_counts.heard);

clear cleanupObj; %#ok<CLCLR>
end

function verifyFile(testCase, filepath)
info = dir(filepath);
testCase.verifyTrue(~isempty(info) && info.bytes > 0);
end
