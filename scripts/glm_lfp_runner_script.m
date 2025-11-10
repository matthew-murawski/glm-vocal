%% Ultra-Compact GLM-LFP Runner

% === EDIT THESE ===



lfp_file = fullfile(P.raw_data_root, 'LFP_Data/M93A/Continuous/M93A_S177_continuousLFP_fixed.mat');
heard = fullfile(P.repo_root, 'data/Label Files/S177/M93A_S177_heard.txt');
produced = fullfile(P.repo_root, 'data/Label Files/S177/M93A_S177_produced.txt');

% build a unique output directory under results/glm with session + timestamp
baseOut = fullfile(P.github_path, 'glm-vocal', 'results', 'glm');
if ~exist(baseOut, 'dir'); mkdir(baseOut); end
tok = regexp(heard, 'S(\d+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(lfp_file, 'S(\d+)', 'tokens', 'once');
end
if ~isempty(tok)
    session_slug = ['s', tok{1}];
else
    session_slug = 's000';
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(baseOut, sprintf('%s_%s', lower(session_slug), timestamp));
% ==================

run_fit_lfp_multichannel('config/defaults.json', lfp_file, heard, produced, outdir);
