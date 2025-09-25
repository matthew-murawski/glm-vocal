function tests = test_build_history_block
% section registration
% we expose the local tests so the matlab runner can pick them up.
tests = functiontests(localfunctions);
end

function stim = makeStim(nT, dt)
% section stimulus helper
% we prepare a minimal stimulus struct that mirrors the production timeline metadata.
stim = struct();
stim.dt = dt;
stim.t = (0:dt:dt * (nT - 1))';
stim.mask = struct('good', true(nT, 1));
end

function testHistoryColumnsAreShiftedSpikes(testCase)
% section causal shift behaviour
% we confirm each column corresponds to progressively deeper spike history without including the current bin.
stim = makeStim(6, 0.1);
sps = [0; 2; 1; 0; 3; 1];
[Xblk, info] = build_history_block(sps, stim, [0.1, 0.3]);

expectedLags = (1:3)';
testCase.verifyEqual(info.mode, 'history');
testCase.verifyEqual(info.lag_bins, expectedLags);
testCase.verifyEqual(info.lag_times_s, expectedLags * stim.dt, 'AbsTol', 0);

fullX = full(Xblk);
testCase.verifyEqual(fullX(:, 1), [0; sps(1:end-1)]);
testCase.verifyEqual(fullX(:, 2), [0; 0; sps(1:end-2)]);
testCase.verifyEqual(fullX(:, 3), [0; 0; 0; sps(1:end-3)]);
end

function testWindowMustExcludeZeroLag(testCase)
% section input validation
% we expect an informative error when the caller specifies a window that tries to include the current bin.
stim = makeStim(5, 0.05);
sps = (1:5)';

verifyError(testCase, @() build_history_block(sps, stim, [0, 0.1]), ...
    'build_history_block:HistoryWindow');
end
