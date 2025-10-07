function demo_s177()
% section overview
% run the single-neuron glm pipeline on all of s177.

[rootDir, ~] = resolve_paths();
dataDir = '/Users/matt/Documents/GitHub/vocalization/data/Label Files/S177';

cfgPath = fullfile(rootDir, 'config', 'defaults.json');
spikePath = fullfile(dataDir, 'M93A_S177_spike_times_ch28.mat');
heardPath = fullfile(dataDir, 'M93A_S177_heard.txt');
producedPath = fullfile(dataDir, 'M93A_S177_produced.txt');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile(rootDir, 'results', ['S177_' timestamp]);

cfg = jsondecode(fileread(cfgPath));
if ~isfield(cfg, 'exclude_predictors') || isempty(cfg.exclude_predictors)
    cfg.exclude_predictors = {'states', 'spike_history'};
end
cfg.produced_split_mode = 'call_type';

fprintf('Running S177 → output: %s\n', outdir);
result = run_fit_single_neuron(cfg, spikePath, heardPath, producedPath, outdir);

fprintf('Finished. Inspect %s for artifacts.\n', outdir);

% section produced-call preview
% reuse the first produced call to show the new prediction plot for the baseline (with history) model.
produced_mask = arrayfun(@(evt) isfield(evt, 'kind') && strcmpi(evt.kind, 'produced'), result.events);
produced_events = result.events(produced_mask);
if ~isempty(produced_events)
    first_produced_call_time = produced_events(1).t_on;
    if ~isempty(first_produced_call_time) && isfinite(first_produced_call_time)
        start_sec = first_produced_call_time - 5;
        duration_sec = 20;

        predicted_series = [];
        if isfield(result.rate, 'mu')
            predicted_series = result.rate.mu;
        end

        if ~isempty(predicted_series)
            plot_predictions(result.stim, result.sps, predicted_series, result.events, start_sec, duration_sec, 'Model Prediction with Spike History');
        end
    end
end
end

function [rootDir, srcDir] = resolve_paths()
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
srcDir = genpath(fullfile(rootDir, 'src'));
addpath(srcDir);
end
