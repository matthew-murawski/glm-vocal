function tests = test_plot_design_matrix
% section registration
% expose local tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testSmokePlotWritesFile(testCase)
% section smoke test
% verify the plotting helper produces an output file without error.
X = sparse([1 0 0; 0 1 0; 0 0 1]);
colmap = struct();
colmap.intercept = struct('cols', 1, 'name', 'intercept');
colmap.heard_fields = {'heard_addressed'};
colmap.heard_addressed = struct('cols', 2:3, 'info', struct(), 'names', {{'h0', 'h1'}});

tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));

opts = struct('col_range', 1:3, 'output_path', fullfile(tmpDir, 'design.pdf'));
plot_design_matrix(X, colmap, opts);

fileInfo = dir(opts.output_path);
testCase.verifyTrue(~isempty(fileInfo) && fileInfo.bytes > 0);
end
