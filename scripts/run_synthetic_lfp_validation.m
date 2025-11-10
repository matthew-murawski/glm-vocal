%RUN_SYNTHETIC_LFP_VALIDATION Click-and-run validation of LFP GLM pipeline.
%   Generates synthetic LFP data with known ground-truth relationships,
%   fits the GLM, and displays results for manual inspection.
%
%   This script validates that the GLM pipeline can recover known event-locked
%   responses by fitting to synthetic data where the ground truth is controlled.

%% Parameters - modify these as needed
params = struct();
params.session_duration = 600;           % seconds
params.lfp_fs = 1000;                    % Hz
params.n_channels = 4;
params.perceived_to_produced_ratio = 5.0;
params.mean_intercall_interval = 3.0;    % seconds
params.call_duration_mean = 0.15;        % seconds
params.call_duration_std = 0.05;
params.bout_probability = 0.7;
params.bout_length_mean = 3;
params.addressed_window = 10.0;          % seconds
params.gamma_center_freq = 110;          % Hz
params.gamma_bandwidth = 40;             % Hz (so 70-150 Hz band)
params.produced_response_amplitude = 2.0;
params.addressed_response_amplitude = 1.0;
params.overheard_response_amplitude = 0.3;
params.perceived_latency = 0.3;          % seconds
params.produced_pre_duration = 0.5;      % seconds
params.noise_amplitude = 0.5;
params.snr_db = 10;

% Output directory
outdir = fullfile(P.github_path, 'glm-vocal/results_synthetic', ['synth_' datestr(now, 'yyyymmdd_HHMMSS')]);

% Config path
cfgPath = fullfile(P.github_path, 'glm-vocal/config', 'defaults.json');

%% Generate synthetic data
fprintf('\n========================================\n');
fprintf('SYNTHETIC LFP DATA GENERATION\n');
fprintf('========================================\n\n');

fprintf('Generating synthetic LFP data...\n');
[lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);
fprintf('\nSynthetic data saved to: %s\n', outdir);

%% Fit GLM to synthetic data
fprintf('\n========================================\n');
fprintf('GLM FITTING\n');
fprintf('========================================\n\n');

fprintf('Fitting GLM to synthetic data...\n');
    results = run_fit_lfp_multichannel(cfgPath, lfpFile, heardFile, producedFile, outdir, P);

%% Display validation summary
fprintf('\n========================================\n');
fprintf('VALIDATION SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Session parameters:\n');
fprintf('  Duration: %.1f s\n', params.session_duration);
fprintf('  Sampling rate: %d Hz\n', params.lfp_fs);
fprintf('  Number of channels: %d\n', params.n_channels);
fprintf('  Gamma band: %.0f-%.0f Hz\n', ...
        params.gamma_center_freq - params.gamma_bandwidth, ...
        params.gamma_center_freq + params.gamma_bandwidth);

fprintf('\nEvent counts:\n');
fprintf('  Perceived calls: %d\n', count_events(heardFile));
fprintf('  Produced calls: %d\n', count_events(producedFile));
fprintf('  Perceived/Produced ratio: %.1f\n', params.perceived_to_produced_ratio);

fprintf('\nFit quality (R² per channel):\n');
r2_values = zeros(length(results.results), 1);
for ch = 1:length(results.results)
    r2_values(ch) = results.results(ch).metrics.r2;
    fprintf('  Channel %d: R² = %.3f', ch, r2_values(ch));
    if r2_values(ch) > 0.3
        fprintf(' [GOOD]\n');
    elseif r2_values(ch) > 0.1
        fprintf(' [FAIR]\n');
    else
        fprintf(' [POOR]\n');
    end
end
fprintf('  Mean R²: %.3f\n', mean(r2_values));

fprintf('\nGround truth expectations:\n');
fprintf('  - Produced kernels should show pre-call activity (before t=0)\n');
fprintf('  - Perceived kernels should peak around t=%.2f s\n', params.perceived_latency);
fprintf('  - Addressed responses should be stronger than overheard\n');
fprintf('  - R² values should be > 0.3 with SNR = %.0f dB\n', params.snr_db);

fprintf('\nValidation outputs:\n');
fprintf('  Plots directory: %s/plots/\n', outdir);
fprintf('    - R² by channel\n');
fprintf('    - LFP traces (actual vs predicted)\n');
fprintf('    - Kernel shapes (lfp_kernels_summary_ch*.pdf)\n');
fprintf('  Results file: %s/fit_results_lfp.mat\n', outdir);

fprintf('\n========================================\n');
fprintf('VALIDATION COMPLETE\n');
fprintf('========================================\n\n');

fprintf('Next steps:\n');
fprintf('  1. Inspect plots in the output directory\n');
fprintf('  2. Verify fitted kernels match ground truth expectations\n');
fprintf('  3. If R² is low, try adjusting SNR or response amplitudes\n');
fprintf('  4. Once validated, apply pipeline to real experimental data\n\n');

%% Helper function
function n = count_events(labelFile)
%COUNT_EVENTS Quick count of events in an Audacity TXT file.
if isempty(labelFile) || ~exist(labelFile, 'file')
    n = 0;
    return;
end

fid = fopen(labelFile, 'r');
if fid < 0
    n = 0;
    return;
end

n = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line)) && ~startsWith(strtrim(line), '#')
        n = n + 1;
    end
end
fclose(fid);
end

