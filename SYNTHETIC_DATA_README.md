# Synthetic LFP Data Generation for GLM Validation

## Overview

This validation system generates synthetic LFP (Local Field Potential) data with known ground-truth event-locked responses, allowing you to verify that the GLM pipeline correctly recovers the relationships you bake into the data.

## Quick Start

1. Open MATLAB and navigate to the repository root directory
2. Run the validation script:
   ```matlab
   run_synthetic_lfp_validation
   ```
3. Check the output in `results_synthetic/synth_YYYYMMDD_HHMMSS/`

## Files Created

### Core Functions

- **`src/io/generate_synthetic_lfp_data.m`** - Main generator function that creates:
  - Synthetic LFP traces with event-locked gamma responses
  - Audacity-format event label files (perceived and produced calls)
  - Realistic conversational structure with bouts

- **`scripts/run_synthetic_lfp_validation.m`** - Click-and-run validation script that:
  - Generates synthetic data with configurable parameters
  - Fits the GLM pipeline to synthetic data
  - Displays validation summary and fit quality metrics

## How It Works

### Event Generation

The generator creates a realistic conversational structure:

1. **Conversational bouts**: Produced calls cluster together with `bout_probability`
2. **Perceived/produced ratio**: Controls the balance (default 5:1)
3. **Addressed vs overheard**: Perceived calls within `addressed_window` (10s) before a produced call are marked as "addressed"

### LFP Signal Generation

For each event type, gamma-band responses are generated:

1. **Produced calls**:
   - Pre-call activity starts `produced_pre_duration` (0.5s) before onset
   - Strong gamma burst (amplitude: 2.0) that spans pre-call through call duration
   - Peaks around call onset

2. **Perceived calls (addressed)**:
   - Response starts at `perceived_latency` (0.3s) after call onset
   - Moderate gamma burst (amplitude: 1.0)
   - Duration ~250ms

3. **Perceived calls (overheard)**:
   - Same timing as addressed
   - Much weaker response (amplitude: 0.3)

### Signal Composition

- **Baseline**: Pink noise (1/f spectrum) for realistic baseline
- **Gamma bursts**: Modulated oscillations at `gamma_center_freq` (100 Hz) ± `gamma_bandwidth` (40 Hz)
- **Envelope**: Gaussian envelope with smooth onset/offset
- **Filtering**: Bandpass filtered to high-gamma range (80-120 Hz by default)
- **SNR**: Adjustable signal-to-noise ratio (default: 10 dB)

## Configuration

### Key Parameters (in `run_synthetic_lfp_validation.m`)

```matlab
params.session_duration = 600;              % Recording length (seconds)
params.lfp_fs = 1000;                       % Sampling rate (Hz)
params.n_channels = 4;                      % Number of LFP channels
params.perceived_to_produced_ratio = 5.0;   % Ratio of heard to produced calls
params.mean_intercall_interval = 3.0;       % Average time between calls (s)
params.addressed_window = 10.0;             % Window for addressed classification (s)

% Response timing
params.perceived_latency = 0.3;             % Delay for perceived response (s)
params.produced_pre_duration = 0.5;         % Pre-call activity duration (s)

% Response amplitudes
params.produced_response_amplitude = 2.0;   % Produced call response strength
params.addressed_response_amplitude = 1.0;  % Addressed perceived strength
params.overheard_response_amplitude = 0.3;  % Overheard perceived strength

% Signal properties
params.gamma_center_freq = 100;             % Center frequency (Hz)
params.gamma_bandwidth = 40;                % Bandwidth (Hz)
params.noise_amplitude = 0.5;               % Baseline noise level
params.snr_db = 10;                         % Signal-to-noise ratio (dB)
```

## Validation Workflow

1. **Generate synthetic data** with known ground truth
2. **Fit GLM** using existing pipeline (`run_fit_lfp_multichannel`)
3. **Inspect results**:
   - Check R² values (should be > 0.3 with SNR=10dB)
   - Verify kernel timing matches expectations
   - Compare addressed vs overheard response magnitudes

### Expected Results

When validation is successful:

- **Produced kernels**: Show pre-call activity (weights before t=0)
- **Perceived kernels**: Peak around t=0.3s (configurable latency)
- **Addressed vs overheard**: Addressed responses should be ~3× stronger (1.0 vs 0.3)
- **R² values**: Should be reasonably high (>0.3) given the SNR

## Output Files

After running `run_synthetic_lfp_validation`, check:

```
results_synthetic/synth_YYYYMMDD_HHMMSS/
├── synthetic_lfp.mat              # LFP data (n_samples × n_channels)
├── synthetic_heard.txt            # Perceived call labels (Audacity format)
├── synthetic_produced.txt         # Produced call labels (Audacity format)
├── fit_results_lfp.mat           # GLM fitting results
└── plots/
    ├── lfp_r2_by_channel.pdf     # R² across channels
    └── lfp_traces_ch*.pdf        # Actual vs predicted traces
```

## Troubleshooting

### Low R² values

- Increase `snr_db` (more signal relative to noise)
- Increase response amplitudes
- Increase `session_duration` or decrease `mean_intercall_interval` (more events)

### Kernels don't match expectations

- Check event files were generated correctly (count events)
- Verify `perceived_latency` and `produced_pre_duration` match your expectations
- Ensure `gamma_center_freq` is in a range your LFP analysis captures

### Pipeline errors

- Verify `config/defaults.json` exists and is valid
- Check that `src/` directory is in MATLAB path
- Ensure all GLM pipeline functions are available

## Advanced Usage

### Custom Parameter Sets

Create multiple validation runs with different parameters:

```matlab
% Test different SNR levels
for snr = [5, 10, 15, 20]
    params.snr_db = snr;
    outdir = sprintf('results_synthetic/snr_%d', snr);
    [lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);
    results = run_fit_lfp_multichannel(cfgPath, lfpFile, heardFile, producedFile, outdir);
end
```

### Using Generator Standalone

Generate data without fitting:

```matlab
addpath('src/io');
params = struct();
params.session_duration = 300;
params.n_channels = 8;
outdir = 'my_synthetic_data';
[lfpFile, heardFile, producedFile] = generate_synthetic_lfp_data(outdir, params);
```

## Implementation Notes

### Design Decisions

1. **Gaussian envelope**: Smooth onset/offset for gamma bursts (could also use raised cosine or Hanning)
2. **Pink noise**: 1/f spectrum baseline (approximated by filtered white noise)
3. **Channel correlation**: Channels have same event responses but independent noise
4. **Frequency jitter**: ±5 Hz per burst for realism

### Helper Functions

The generator includes modular helper functions:
- `generate_event_times()` - Creates conversational structure
- `classify_addressed_calls()` - Marks addressed vs overheard
- `generate_synthetic_lfp()` - Creates LFP traces
- `generate_gamma_burst()` - Single gamma response
- `generate_baseline_noise()` - Pink noise baseline
- `write_audacity_labels()` - Exports to Audacity format

## Next Steps

1. Run validation to confirm pipeline works on synthetic data
2. If successful, apply to real experimental data
3. Adjust synthetic data parameters to match properties of your real recordings
4. Use as regression test when modifying the GLM pipeline

## Questions?

If anything is unclear or you encounter issues:
- Check that all parameters are sensible for your use case
- Verify file paths and directory permissions
- Inspect intermediate outputs (event files, LFP data) manually
- Compare synthetic data properties to real recordings
