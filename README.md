# GLM Vocal

Repo for fitting a generalized linear model (GLM) to natural marmoset conversation data. The pipeline ingests spike times alongside produced and perceived vocalizations, rasterizes them onto a common grid, constructs lagged design matrices, and fits a Poisson GLM with smoothness penalties to recover stimulus kernels, spike-history effects, and conversational context weights.

## What This Repository Provides
- End-to-end MATLAB pipeline (`scripts/run_fit_single_neuron.m`) that loads spikes + labels, builds regressors, fits the GLM, and saves metrics, plots, and artifacts.
- Preprocessing utilities for timebase construction, spike binning, conversational state detection, and regressor stream generation.
- Feature builders that assemble sparse design matrices with causal/symmetric kernels and smoothness penalties.
- Model routines for penalized MAP fitting, blocked cross-validation, rate prediction, and parameter unpacking.
- Evaluation helpers and plotting scripts for quick QC (kernels, rate vs. spikes, CV curves, PSTHs, design matrix snapshots).
- Unit and end-to-end tests runnable from MATLAB (or via the provided shell wrapper) to keep the stack stable.

## Repository Layout
```
glm-vocal/
├─ README.md                     # project overview (this file)
├─ docs/                         # high-level blueprint and spec
├─ config/                       # JSON configuration defaults
├─ demos/                        # small demo dataset for S178 session snippet
├─ scripts/                      # entry points (fits, demos, exports, tests)
├─ src/                          # MATLAB source grouped by subsystem
├─ tests/                        # unit and e2e MATLAB tests
└─ results/                      # example fit artifacts (outputs land here by default)
```

See `docs/SPEC.md` for the formal MVP requirements and `docs/BLUEPRINT.md` for the phased implementation plan.

## Requirements
- MATLAB R2021b or newer (tested with recent desktop MATLAB builds).
- Optimization Toolbox (the code relies on `fminunc`).
- macOS paths are assumed in a few helper scripts; adjust to suit your environment if needed.

## Quick Start
### 1. Clone and open in MATLAB
Add the repository to your MATLAB path, or run scripts directly with `matlab -batch ...`.

### 2. Run the demo dataset (optional)
A 60 s snippet from session S178 is bundled under `demos/data`. From MATLAB:
```matlab
scripts.demo_s178_snippet
```
This resolves the config, demo spike file, and labels, then writes artifacts to `results/demo_s178_<timestamp>/`.

### 3. Fit your own session
`run_fit_single_neuron` expects spike times (MAT file) plus optional heard/produced labels (MAT or Audacity TXT):
```matlab
out = scripts.run_fit_single_neuron( ...
    "config/defaults.json", ...      % configuration file (fallbacks to defaults if empty)
    "data/sessionA/neuron001.mat", ... % spike MAT with spike_times, neuron_id, session_id
    "data/sessionA/heard.txt", ...     % perceived-call labels (optional)
    "data/sessionA/produced.txt", ...  % produced-call labels (optional)
    "results/neuron001" ...            % output directory (auto-created)
);
```
Arguments after the spike file are optional; pass `[]` to fall back to defaults. The function returns a struct with configuration, fit results, metrics, and output paths.

### Input expectations
- **Spikes (`load_spikes`)**: MAT file containing `spike_times` (double vector, seconds), `neuron_id`, and `session_id`.
- **Labels (`load_labels`)**: Either MAT struct arrays with fields `kind`, `onset`, `offset`, `label` or Audacity TXT exports (`start\tstop\toptional_label`). `kind` must be `produced` or `perceived`.
- **Config (`config/defaults.json`)**: Controls bin size (`dt`), kernel windows, cross-validation grid, state rules, optimizer tolerances, and RNG seed.

## Outputs
`run_fit_single_neuron` writes a timestamped subdirectory (default `results/`) containing:
- `fit_results.mat` with the configuration, preprocessed streams, design matrix, smoothness operators, selected λ, fitted weights, metrics, and QC summary.
- Plots: design matrix preview, kernel traces, firing-rate vs spikes, CV curve, PSTHs (stored under `results/<run>/plots/`).
- JSON/text summary files from `qc_session_summary` describing the fit quality and artifact locations.

## Testing
Run the full MATLAB test suite (unit + e2e) from the repo root:
```
./scripts/run_matlab_tests.sh
```
This wrapper adds `src/` and `tests/` to the MATLAB path and executes all suites with `runtests`.

## Development Notes
- The pipeline is organized to encourage sparse operations and deterministic workflows; see `docs/BLUEPRINT.md` for phased milestones.
- All MATLAB functions follow the comment style described in the agent instructions (section headings, informal tone).
- Additional helper scripts (`scripts/export_state_labels_to_audacity.m`, `scripts/demo_s178.m`) show how to adapt the main entry point to different datasets or export formats.

## Questions & Further Reading
- **Detailed specification**: `docs/SPEC.md`
- **Implementation plan**: `docs/BLUEPRINT.md`
