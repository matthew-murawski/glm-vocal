%% run_eta_sanity_check
% event-triggered average sanity check for high-gamma power

clr; clc;

%% User Parameters
animalID = 'M93A';
sessionNumber = 177;

dataPath = fullfile(P.repo_root, 'output/glm/M93A_S177_HG_for_GLM.mat');
etaWindow_s = [-2.0, 2.0];
target_dt = 0.1;

ChanGeo = [ 3,  7, 11, 15, 17, 21, 25, 29; 
            1,  5,  9, 13, 19, 23, 27, 31; 
            4,  8, 12, 16, 18, 22, 26, 30; 
            2,  6, 10, 14, 20, 24, 28, 32];

%% Execute Analysis
fprintf('Starting ETA sanity check for %s | session %d\n', animalID, sessionNumber);
orchestrate_eta_sanity_check(animalID, sessionNumber, dataPath, ChanGeo, etaWindow_s, target_dt);
fprintf('ETA sanity check complete.\n');