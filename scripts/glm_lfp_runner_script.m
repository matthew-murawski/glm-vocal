%% Ultra-Compact GLM-LFP Runner

% === EDIT THESE ===
lfp_file = 'M93A_S177_continuousLFP_fixed.mat';
heard = 'M93A_S177_heard.txt';
produced = 'M93A_S177_produced.txt';
outdir = 'results/M93A_S177';
% ==================

run_fit_lfp_multichannel('config/defaults.json', lfp_file, heard, produced, outdir);