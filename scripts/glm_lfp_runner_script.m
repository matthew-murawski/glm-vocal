%% Ultra-Compact GLM-LFP Runner

% === EDIT THESE ===



lfp_file = fullfile(P.raw_data_root, 'LFP_Data/M93A/Continuous/M93A_S177_continuousLFP_fixed.mat');
heard = fullfile(P.repo_root, 'data/Label Files/S177/M93A_S177_heard.txt');
produced = fullfile(P.repo_root, 'data/Label Files/S177/M93A_S177_produced.txt');
outdir = fullfile(P.github_path, 'glm-vocal/results/M93A_S177');
% ==================

run_fit_lfp_multichannel('config/defaults.json', lfp_file, heard, produced, outdir);