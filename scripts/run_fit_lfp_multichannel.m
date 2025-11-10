%% run_fit_lfp_multichannel
% one-click runner for multichannel LFP GLM fitting

clr; clc;

%% User Parameters
monkey = 'M93A';
session = 177;

%% Execute Analysis
fprintf('Starting LFP GLM fit for %s | session %d\n', monkey, session);

% build absolute paths from environment struct P
cfgPath = fullfile(P.github_path, 'glm-vocal', 'config', 'defaults.json');
% lfpPath = fullfile(P.raw_data_root, 'LFP_MAT', monkey, 'Continuous', ...
%                   sprintf('%s_S%d_continuousLFP.mat', monkey, session));

lfpPath = fullfile(P.repo_root, 'output/glm/M93A_S177_HG_for_GLM.mat');

labelDir = fullfile(P.repo_root, 'data', 'Label Files', sprintf('S%d', session));
heardPath = fullfile(labelDir, sprintf('%s_S%d_heard.txt', monkey, session));
producedPath = fullfile(labelDir, sprintf('%s_S%d_produced.txt', monkey, session));
% make a unique output folder under results/glm using session and timestamp
baseOut = fullfile(P.github_path, 'glm-vocal', 'results', 'glm');
if ~exist(baseOut, 'dir'); mkdir(baseOut); end
session_slug = sprintf('s%d', session);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(baseOut, sprintf('%s_%s', lower(session_slug), timestamp));
if ~exist(outdir, 'dir'); mkdir(outdir); end

% run orchestrator
results = orchestrate_fit_lfp_multichannel(cfgPath, lfpPath, heardPath, producedPath, outdir, P);

fprintf('Results generated at %s. Done.\n', outdir);
