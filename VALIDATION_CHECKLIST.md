# Synthetic LFP Validation Checklist

Use this checklist when running the synthetic data validation for the first time.

## Pre-Flight Checks

- [ ] MATLAB is installed and working
- [ ] Repository root is the current working directory
- [ ] `config/defaults.json` exists
- [ ] `src/` directory is in MATLAB path (or will be added automatically)

## Step 1: Generate Synthetic Data Only

Test the generator independently before running the full pipeline:

```matlab
% Navigate to repository root
cd('/path/to/glm-vocal')

% Add source to path
addpath('src/io')

% Create test parameters (short session for quick test)
params = struct();
params.session_duration = 60;     % 1 minute
params.n_channels = 2;            % Just 2 channels
params.lfp_fs = 1000;

% Generate data
outdir = 'test_synthetic_quick';
[lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);
```

### Verify Outputs

- [ ] Directory created: `test_synthetic_quick/`
- [ ] LFP file exists: `test_synthetic_quick/synthetic_lfp.mat`
- [ ] Heard file exists: `test_synthetic_quick/synthetic_heard.txt`
- [ ] Produced file exists: `test_synthetic_quick/synthetic_produced.txt`

### Inspect LFP File

```matlab
% Load and check LFP data
data = load(lfpFile);
whos('-file', lfpFile)

% Verify structure
assert(isfield(data, 'lfp_data'), 'Missing lfp_data');
assert(isfield(data, 'lfp_fs'), 'Missing lfp_fs');
assert(isfield(data, 'lfp_t0'), 'Missing lfp_t0');
assert(isfield(data, 'channel_ids'), 'Missing channel_ids');

% Check dimensions
fprintf('LFP data size: %d samples × %d channels\n', size(data.lfp_data));
fprintf('Expected: ~%d samples × %d channels\n', 60*1000, 2);

% Check values are reasonable
fprintf('LFP data range: [%.2f, %.2f]\n', min(data.lfp_data(:)), max(data.lfp_data(:)));
fprintf('LFP data mean: %.4f (should be near 0)\n', mean(data.lfp_data(:)));
```

Expected output:
- Size: ~60000 × 2
- Range: reasonable values (not all zeros, no Inf/NaN)
- Mean: close to 0

### Inspect Event Files

```matlab
% Read heard events
heardEvents = load_labels(heardFile, 'perceived');
fprintf('Perceived calls: %d\n', numel(heardEvents));
fprintf('First 5 events:\n');
for ii = 1:min(5, numel(heardEvents))
    fprintf('  t=%.3f-%.3f: %s\n', heardEvents(ii).t_on, ...
            heardEvents(ii).t_off, heardEvents(ii).label);
end

% Read produced events
producedEvents = load_labels(producedFile, 'produced');
fprintf('\nProduced calls: %d\n', numel(producedEvents));
fprintf('First 5 events:\n');
for ii = 1:min(5, numel(producedEvents))
    fprintf('  t=%.3f-%.3f\n', producedEvents(ii).t_on, producedEvents(ii).t_off);
end

% Check ratio
ratio = numel(heardEvents) / numel(producedEvents);
fprintf('\nPerceived/Produced ratio: %.2f (expected: 5.0)\n', ratio);
```

Expected output:
- ~10-20 perceived calls
- ~2-4 produced calls
- Ratio ~5:1 (may vary due to randomness)
- Labels: "addressed" or "overheard"

## Step 2: Run Full Validation Pipeline

Once the generator works, run the full validation:

```matlab
% Run full validation (this will take longer)
run_synthetic_lfp_validation
```

### Check Console Output

Expected output sections:

1. **Data Generation**
   - [ ] "Generating conversational event structure..."
   - [ ] "Generating synthetic LFP data..."
   - [ ] "Saving synthetic data files..."

2. **GLM Fitting**
   - [ ] "Fitting GLM to X channels..."
   - [ ] R² values printed for each channel

3. **Validation Summary**
   - [ ] Event counts displayed
   - [ ] R² values > 0.3 marked as [GOOD]
   - [ ] Ground truth expectations listed

### Check Output Files

- [ ] Results directory created: `results_synthetic/synth_YYYYMMDD_HHMMSS/`
- [ ] Plots directory exists: `results_synthetic/synth_YYYYMMDD_HHMMSS/plots/`
- [ ] R² plot exists: `plots/lfp_r2_by_channel.pdf`
- [ ] Trace plots exist: `plots/lfp_traces_ch*.pdf`
- [ ] Results MAT file exists: `fit_results_lfp.mat`

### Inspect Plots

1. **R² by Channel** (`lfp_r2_by_channel.pdf`)
   - [ ] All channels have R² > 0.1
   - [ ] Most channels have R² > 0.3
   - [ ] Values are consistent across channels

2. **LFP Traces** (`lfp_traces_ch*.pdf`)
   - [ ] Predicted trace (orange/red) follows actual trace (blue)
   - [ ] Peaks align with event times
   - [ ] Not just fitting noise

### Inspect Kernels

```matlab
% Load results
load('results_synthetic/synth_*/fit_results_lfp.mat');

% Plot produced kernel for channel 1
ch = 1;
if isfield(results(ch).kernels, 'produced_any')
    kernel = results(ch).kernels.produced_any;
    figure;
    plot(kernel.t, kernel.h);
    xlabel('Time (s)');
    ylabel('Kernel weight');
    title('Produced Call Kernel');
    grid on;

    % Check for pre-call activity
    pre_call_idx = kernel.t < 0;
    if any(abs(kernel.h(pre_call_idx)) > 0.1)
        fprintf('✓ Produced kernel shows pre-call activity\n');
    else
        fprintf('✗ WARNING: Produced kernel missing pre-call activity\n');
    end
end

% Plot perceived kernel
if isfield(results(ch).kernels, 'heard_any')
    kernel = results(ch).kernels.heard_any;
    figure;
    plot(kernel.t, kernel.h);
    xlabel('Time (s)');
    ylabel('Kernel weight');
    title('Perceived Call Kernel');
    grid on;

    % Check for post-call peak around 0.3s
    [max_val, max_idx] = max(abs(kernel.h));
    peak_time = kernel.t(max_idx);
    fprintf('Perceived kernel peak at t=%.3f s (expected: ~0.3 s)\n', peak_time);

    if abs(peak_time - 0.3) < 0.1
        fprintf('✓ Perceived kernel peak timing is correct\n');
    else
        fprintf('✗ WARNING: Perceived kernel peak at wrong time\n');
    end
end
```

## Step 3: Validation Criteria

### Success Criteria

- [x] **Data generation completes without errors**
- [ ] **Event files are valid Audacity format**
- [ ] **LFP data has correct dimensions**
- [ ] **GLM fitting completes without errors**
- [ ] **Mean R² > 0.3 across channels**
- [ ] **Produced kernels show pre-call activity**
- [ ] **Perceived kernels peak around t=0.3s**
- [ ] **Addressed responses stronger than overheard** (if kernel splitting is enabled)

### If Validation Fails

#### Low R² values (< 0.2)

Possible causes:
- SNR too low → increase `params.snr_db`
- Response amplitudes too weak → increase amplitude parameters
- Not enough events → increase `session_duration` or decrease `mean_intercall_interval`
- Baseline noise too strong → decrease `noise_amplitude`

Try:
```matlab
params.snr_db = 20;  % increase from 10
params.produced_response_amplitude = 3.0;  % increase from 2.0
```

#### Kernel timing is wrong

Possible causes:
- Config `dt` doesn't match expectations
- Latency parameters set incorrectly
- Sampling rate mismatch

Check:
```matlab
% Verify config dt
cfg = jsondecode(fileread('config/defaults.json'));
fprintf('Config dt: %.4f s\n', cfg.dt);
```

#### Pipeline errors

Common issues:
- Missing `config/defaults.json` → verify it exists
- Path issues → ensure `src/` is in MATLAB path
- Missing dependencies → check all functions load correctly

## Step 4: Parameter Sweep (Optional)

Test robustness by varying parameters:

```matlab
% Test different SNR levels
for snr_db = [5, 10, 15, 20]
    params = struct();
    params.session_duration = 300;
    params.snr_db = snr_db;

    outdir = sprintf('results_synthetic/snr_%d', snr_db);
    [lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);
    results = run_fit_lfp_multichannel('config/defaults.json', ...
                                        lfpFile, heardFile, producedFile, outdir);

    % Compute mean R²
    r2_vals = arrayfun(@(r) r.metrics.r2, results.results);
    fprintf('SNR=%d dB: Mean R²=%.3f\n', snr_db, mean(r2_vals));
end
```

Expected: R² should increase with SNR

## Step 5: Compare to Real Data

Once validation succeeds:

1. Run GLM on real experimental data
2. Compare:
   - R² values (real vs synthetic)
   - Kernel shapes and timing
   - Response amplitudes
3. Adjust synthetic parameters to match real data properties
4. Use as regression test for pipeline changes

## Troubleshooting

### Error: "missing required field X"

Check that the MAT file or event files have correct structure. Run Step 1 checks.

### Error: "unable to open file"

Check file paths and permissions. Ensure output directory is writable.

### Warning: "LFP data was transposed"

Normal if channels > samples. Verify dimensions are correct after loading.

### Plots look noisy / no clear signal

- Increase response amplitudes
- Increase SNR
- Check that gamma_center_freq matches your filter settings

### Kernels are flat / no structure

- Check that events were generated correctly (non-empty)
- Verify response amplitudes are non-zero
- Ensure event times fall within session duration

---

## Quick Command Reference

```matlab
% Generate data only
[lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);

% Run full validation
run_synthetic_lfp_validation

% Load and inspect results
load('results_synthetic/synth_*/fit_results_lfp.mat');
load('results_synthetic/synth_*/synthetic_lfp.mat');

% Count events
load_labels('results_synthetic/synth_*/synthetic_heard.txt', 'perceived');
load_labels('results_synthetic/synth_*/synthetic_produced.txt', 'produced');
```
