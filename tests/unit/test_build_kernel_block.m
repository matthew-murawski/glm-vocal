function tests = test_build_kernel_block
% section registration
% we expose the test cases to the MATLAB runner so the harness can execute them.
tests = functiontests(localfunctions);
end

function stim = makeStim(nT, dt)
% section stimulus helper
% we construct a minimal stimulus struct that mimics the timeline metadata used in production code.
stim = struct();
stim.dt = dt;
stim.t = (0:dt:dt * (nT - 1))';
stim.mask = struct('good', true(nT, 1));
end

function testColumnsAreShiftedCopies(testCase)
% section core behaviour
% we confirm each column equals the regressor shifted by the corresponding lag with leading zeros applied.
stim = makeStim(6, 0.1);
z = [0; 1; 0; 1; 1; 0];
[Xblk, info] = build_kernel_block(z, stim, [0, 0.3]);

expectedLags = (0:3)';
testCase.verifyEqual(size(Xblk), [numel(z), numel(expectedLags)]);
testCase.verifyEqual(info.lag_bins, expectedLags);
testCase.verifyEqual(info.lag_times_s, expectedLags * stim.dt, 'AbsTol', 0);

for k = 0:numel(expectedLags) - 1
    col = full(Xblk(:, k + 1));
    shifted = zeros(size(z));
    tailLength = numel(z) - k;
    shifted(k + 1:end) = z(1:tailLength);
    testCase.verifyEqual(col, shifted);
end
end

function testZeroPaddingWhenLagExceedsHistory(testCase)
% section extreme lag behaviour
% we request a window that extends beyond the available samples and expect trailing columns to be all zeros.
stim = makeStim(4, 0.1);
z = [1; 0; 1; 0];
[Xblk, info] = build_kernel_block(z, stim, [0, 0.4]);

testCase.verifyEqual(info.lag_bins, (0:4)');

fullX = full(Xblk);
testCase.verifyEqual(fullX(:, 5), zeros(4, 1));
testCase.verifyEqual(fullX(1, 3), 0);
testCase.verifyEqual(fullX(2, 3), 0);
testCase.verifyEqual(fullX(3, 3), z(1));
end

function testSymmetricModeAlignment(testCase)
% section symmetric behaviour
% we ask for a symmetric window with pre and post lags and ensure the block aligns around lag zero as expected.
stim = makeStim(7, 0.1);
z = (1:7)';
[Xblk, info] = build_kernel_block(z, stim, [-0.2, 0.3], 'symmetric');

expectedLags = (-2:3)';
testCase.verifyEqual(info.mode, 'symmetric');
testCase.verifyEqual(info.lag_bins, expectedLags);

zeroIdx = find(info.lag_bins == 0, 1);
testCase.verifyNotEmpty(zeroIdx);
testCase.verifyEqual(full(Xblk(:, zeroIdx)), z);

for c = 1:numel(expectedLags)
    lag = expectedLags(c);
    col = full(Xblk(:, c));
    expected = zeros(size(z));
    for r = 1:numel(z)
        src = r - lag;
        if src >= 1 && src <= numel(z)
            expected(r) = z(src);
        end
    end
    testCase.verifyEqual(col, expected);
end
end
