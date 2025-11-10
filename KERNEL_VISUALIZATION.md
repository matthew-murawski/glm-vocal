# GLM-LFP Kernel Visualizations

This document describes the kernel visualization features added to the GLM-LFP pipeline.

## Overview

The pipeline now automatically generates visualizations of the fitted predictive kernels, showing the temporal dynamics of high gamma band power responses to different vocal events.

## Generated Plots

After running the GLM fitting pipeline, the following kernel visualizations are created:

### 1. Kernel Summary Plots
**Location**: `<outdir>/plots/lfp_kernels_summary_ch<N>.pdf`

Multi-panel figure showing all fitted kernels for a given channel on a single page.

**Features**:
- One subplot per kernel (produced, heard addressed, heard overheard, etc.)
- Color-coded by event type:
  - **Red**: Produced kernels
  - **Blue**: Heard addressed kernels
  - **Green**: Heard overheard kernels
  - **Purple**: History kernels
- Vertical dashed line at t=0 marks event onset
- Channel ID and fit quality (R², correlation) in title
- Time axis in seconds relative to event

**Interpretation**:
- **Produced kernels**: Should show pre-call activity (negative times) ramping up before vocalization onset
- **Heard addressed kernels**: Should show stronger post-stimulus response (~300ms after call)
- **Heard overheard kernels**: Should show weaker post-stimulus response than addressed
- Kernel amplitude indicates contribution to high gamma band power

### 2. Detailed Individual Kernel Plots
**Location**: `<outdir>/plots/channel_<N>/kernels.pdf`

Detailed view of each kernel with confidence intervals (if computed).

**Features**:
- Individual plots for each predictor
- Optional confidence intervals (95% CI)
- Significance indication (dashed lines for non-significant kernels, p ≥ 0.05)
- State weights shown as bar plots
- Intercept shown as horizontal line

## How to Use

### Automatic Generation
Kernel plots are generated automatically when you run:

```matlab
run_fit_lfp_multichannel(cfgPath, lfpFile, heardFile, producedFile, outdir, P)
```

or the validation script:

```matlab
run_synthetic_lfp_validation
```

### Manual Generation
To generate kernel plots for existing fit results:

```matlab
% Load results
load('results/fit_results_lfp.mat', 'results');

% Plot summary for a specific channel
fig = plot_lfp_kernels_summary(results(1).kernels, ...
                               results(1).channel_id, ...
                               results(1).metrics);
saveas(fig, 'my_kernels.pdf');

% Plot detailed kernels
plot_kernels(results(1).kernels, [], 'output_directory');
```

## Expected Kernel Shapes

Based on the ground truth in synthetic data:

### Produced Kernels
- **Time window**: -2 to +3 seconds relative to call onset
- **Expected shape**: Ramp-up starting ~500ms before call (t=-0.5), peaking around call onset (t=0)
- **Biological interpretation**: Motor preparation and execution of vocalization

### Heard Addressed Kernels
- **Time window**: 0 to +2 seconds after call onset
- **Expected shape**: Peak around 300ms post-stimulus
- **Biological interpretation**: Auditory processing of socially-relevant vocalizations

### Heard Overheard Kernels
- **Time window**: 0 to +2 seconds after call onset
- **Expected shape**: Similar timing to addressed but ~3x weaker amplitude
- **Biological interpretation**: Reduced response to non-directed vocalizations

## Customization

### Colors
Edit color scheme in `src/viz/plot_lfp_kernels_summary.m`:

```matlab
color_produced = [0.8, 0.2, 0.2];      % red
color_addressed = [0.2, 0.6, 0.8];     % blue
color_overheard = [0.4, 0.7, 0.4];     % green
```

### Layout
Modify tile layout parameters:

```matlab
n_cols = min(3, n_kernels);  % max 3 columns
n_rows = ceil(n_kernels / n_cols);
```

### Line Width
Adjust kernel line width:

```matlab
plot(ax, times, weights, 'LineWidth', 2.5, 'Color', color);
```

## Troubleshooting

**No kernels plotted**: Check that your config includes the relevant predictors and they weren't excluded via `cfg.exclude_predictors`.

**Kernels look flat**: May indicate insufficient data, poor SNR, or no true relationship. Check R² values and event counts.

**Time axis unexpected**: Verify `cfg.produced_window_s` and `cfg.heard_window_s` match your experimental design.

## References

- Main fitting script: `scripts/run_fit_lfp_multichannel.m`
- Summary plot function: `src/viz/plot_lfp_kernels_summary.m`
- Detailed plot function: `src/plot/plot_kernels.m`
- Kernel extraction: `src/model/unpack_params.m`
