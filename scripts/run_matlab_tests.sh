#!/usr/bin/env bash
set -euo pipefail
# always compute repo root and pass absolute paths to matlab
REPO_ROOT="$(git rev-parse --show-toplevel)"

/Applications/MATLAB_R2025b.app/bin/matlab -batch "reporoot='${REPO_ROOT}'; try, cd(reporoot); catch, end; addpath(genpath(fullfile(reporoot,'src'))); addpath(genpath(fullfile(reporoot,'tests'))); import matlab.unittest.TestSuite; import matlab.unittest.TestRunner; suite = TestSuite.fromFolder(fullfile(reporoot,'tests'),'IncludingSubfolders',true); if isempty(suite), disp('No tests discovered.'); error('NoTestsFound'); end; runner = TestRunner.withTextOutput('Verbosity',3); results = runner.run(suite); disp(results); if any([results.Failed]), error('TestsFailed'); end;"
