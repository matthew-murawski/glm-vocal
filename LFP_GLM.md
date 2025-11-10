# LFP GLM (High‑Gamma) Pipeline

This document explains the end‑to‑end LFP GLM pipeline for high‑gamma power, from inputs (LFP + labels) to fitted weights, predictions, and plots. It focuses on the multichannel runner and the helper functions (“minions”) it orchestrates, describing what each takes in, how it transforms inputs, and what it outputs.

If you’re new to the repo, skim the TL;DR first, then use the function references to dive deeper.

---

## TL;DR — Data Journey

- Labels → load `produced` and `perceived` call intervals with `load_labels` and optionally consolidate twitter bouts.
- Timebase → build uniform time grid `stim.t` at `cfg.dt` via `build_timebase`.
- LFP ingest → `load_lfp` reads the continuous signal; `bin_lfp` bandpasses to high‑gamma, extracts Hilbert envelope, smooths, resamples to `stim.t`, and normalizes.
- Streams & states → `build_streams` rasterizes heard/produced onsets; `compute_states` flags conversational (`convo`) vs spontaneous (`spon`).
- Design matrix → `assemble_design_matrix` concatenates intercept | heard kernels (causal) | produced kernels (symmetric) | state scalars | LFP history block.
- Penalty → `smoothness_penalty` builds a 2nd‑difference operator per kernel block for λ‑regularized smoothness.
- Fit & CV → `crossval_blocked` selects λ (MSE on held‑out folds); `fit_glm_gauss` fits the Gaussian GLM with identity or log link.
- Predict & evaluate → `predict_lfp` gives the channel’s predicted LFP; `metrics_lfp` computes R², corr, RMSE.
- Unpack & plot → `unpack_params` slices weights into named kernels; plots summarize R², traces, and kernel shapes per channel.

---

## Main Runner

- Entry: `scripts/run_fit_lfp_multichannel.m`
- Minimal wrapper: `scripts/glm_lfp_runner_script.m` (edit paths and call the runner)

What it does per channel:
1) loads config and labels, builds the timebase and binary streams/states.
2) ingests the LFP file, extracts high‑gamma power if configured, resamples and normalizes.
3) assembles the design matrix and penalty; cross‑validates λ; fits the GLM.
4) predicts LFP, computes metrics, unpacks kernels, and writes plots + a results MAT.

Inputs
- `cfgPath`: JSON (see `config/defaults.json` for LFP knobs)
- `lfpPath`: MAT with `lfp_data`, `lfp_fs`, `lfp_t0`, optional `channel_ids`
- `heardPath` / `producedPath`: Audacity `.txt` or MAT labels
- `outdir`: output directory

Outputs
- Plots under `<outdir>/plots/` (R² by channel, traces, kernels)
- Results at `<outdir>/fit_results_lfp.mat` (config, LFP meta, design, weights, preds, metrics)

Example usage
- Edit `scripts/glm_lfp_runner_script.m` and run it in MATLAB.

---

## The Minions (Helper Functions)

Below is the ordered list of helpers the runner calls, with inputs, transformations, and outputs.

1) Load labels
- `src/io/load_labels.m`
  - In: path to MAT struct array or Audacity `.txt`; optional default kind.
  - Transform: normalize to `struct('kind','produced|perceived','t_on','t_off','label')`, validate ordering and finiteness.
  - Out: `events` struct array.
- `src/preprocess/consolidate_twitter_bouts.m` (optional)
  - In: events, `boutWindow` seconds
  - Transform: merge rapid twitter syllables into bouts
  - Out: updated events

2) Timebase
- `src/preprocess/build_timebase.m`
  - In: events, dummy spike struct (unused for LFP), `cfg.dt`
  - Transform: set session duration from latest event; build `stim.t = (0:dt:T)'`, `stim.dt`, `stim.mask.good`
  - Out: `stim`

3) LFP ingest and high‑gamma power
- `src/io/load_lfp.m`
  - In: MAT with `lfp_data` (n_samples×n_channels or transposed), `lfp_fs`, `lfp_t0`
  - Transform: validate, orient to (n_samples×n_channels), compute meta fields
  - Out: `lfp` struct (`data`, `fs`, `t0`, `channel_ids`, `n_samples`, `n_channels`, `duration`)
- `src/preprocess/bin_lfp.m`
  - In: `lfp`, `stim`, `cfg`
  - Transform pipeline:
    - Optional bandpass to high‑gamma: `cfg.lfp.bandpass_hz = [low, high]` (defaults `[70, 150]`)
    - Optional high‑gamma envelope: `cfg.lfp.extract_band_power = true` uses Hilbert amplitude; optional smoothing `cfg.lfp.power_smoothing_ms`
    - Linear interpolation onto `stim.t` (out‑of‑range extrapolation → 0)
    - Per‑channel normalization: `cfg.lfp.normalize ∈ {'zscore','minmax','none'}`
  - Out: `lfp_binned` (n_bins×n_channels), aligned to `stim.t`

4) Streams and states
- `src/preprocess/build_streams.m`
  - In: events, `stim`, `cfg`
  - Transform: heard streams `heard_any`, `heard_addressed`, `heard_overheard`; produced streams split by `cfg.produced_split_mode ∈ {'context','call_type'}`; rasterized onset impulses onto `stim.t`.
  - Out: `streams` with `.heard_fields` and `.produced_fields`
- `src/preprocess/compute_states.m`
  - In: events, `stim`, `cfg.state` (e.g., `response_window_s`, `max_seq_gap_s`)
  - Transform: produce conversation intervals using reply‑window chaining; rasterize to `states.convo`; complement to `states.spon`
  - Out: `states`

5) Design matrix assembly
- `src/features/assemble_design_matrix.m`
  - In: `streams`, `states`, one response vector (`lfp_ch`), `cfg`, `stim`
  - Transform (column‑wise concatenation):
    - Intercept: constant 1
    - Heard kernels (causal): basis expansion over `cfg.heard_window_s` with `cfg.heard_basis`
    - Produced kernels (symmetric): basis over `cfg.produced_window_s` with `cfg.produced_basis`
    - States: two scalar columns `convo`, `spon`
    - LFP history block (causal): `lfp_history` over `cfg.lfp.history_window_s`; if `cfg.lfp.history_basis.kind` ≠ `raw`, project raw lags to a raised‑cosine basis
  - Mask: drop rows where `~stim.mask.good`
  - Out: `Xd = struct('X','y','colmap','response_type')`
    - Note: `colmap.lfp_history.cols` indexes the history block; heard/produced blocks include basis info and lag times.

6) Smoothness penalty
- `src/features/smoothness_penalty.m`
  - In: `Xd.colmap`, `cfg` (block enable flags)
  - Transform: build 2nd‑difference operator for each penalized block; concatenate into sparse `D`
  - Out: `D`, `Dmap`
  - Note: by default, heard and produced kernels are penalized; `lfp_history` is not penalized unless you add `cfg.lfp_history = true` (block names must match colmap field names).

7) Cross‑validation over λ
- `src/model/crossval_blocked.m`
  - In: `fitfun`, `Xd`, `D`, `cfg.cv` (`k`, `lambdas`)
  - Transform: contiguous K‑fold split; for each λ, fit on K−1 folds and score held‑out with mean‑squared error (LFP response)
  - Out: `best_lambda`, `cvinfo` (per‑λ per‑fold metrics)

8) Fit the GLM (Gaussian)
- `src/model/fit_glm_gauss.m`
  - In: `Xd`, `D`, `lambda`, `cfg.optimizer`, `cfg.lfp.link ∈ {'identity','log'}`
  - Transform: minimize 0.5·||y − μ||² + λ·||D w||² with μ = Xw (identity) or μ = exp(Xw) (log). Closed‑form ridge is used for identity when possible; otherwise `fminunc`/fallback.
  - Out: `wmap.w` (weights), `wmap.hessian`, `fitinfo`

9) Predict and evaluate
- `src/model/predict_lfp.m`
  - In: `Xd.X`, `wmap.w`, link
  - Out: `lfp_pred`
- `src/eval/metrics_lfp.m`
  - In: `Xd.y`, `lfp_pred`
  - Out: `metrics` (MSE, RMSE, R², corr, explained variance)

10) Unpack and visualize
- `src/model/unpack_params.m`
  - In: `wmap`, `Xd.colmap`, `cfg`, `stim`
  - Out: `kernels` struct with named blocks (intercept, heard_*, produced_*, states, spike_history)
  - Note: the LFP history block is named `lfp_history` in `colmap`, but current unpacker only auto‑extracts `spike_history`. You can still access history weights directly via `Xd.colmap.lfp_history.cols` or extend the unpacker to include it.
- Plots
  - `src/viz/plot_lfp_r2_by_channel.m` — bar chart of R² across channels
  - `src/viz/plot_lfp_traces.m` — actual vs predicted traces per channel
  - `src/viz/plot_lfp_kernels_summary.m` — produced/heard/state kernels per channel
  - `src/plot/plot_kernels.m` — generic multi‑kernel panel (driven by fields present in `kernels`)

---

## Configuration Knobs (High‑Gamma)

Key fields in `config/defaults.json` (override as needed):
- Timebase: `dt` (s) — choose temporal granularity for `stim.t`
- High‑gamma band: `lfp.bandpass_hz` (e.g., `[70,150]`)
- Envelope extraction: `lfp.extract_band_power = true` (Hilbert amplitude)
- Envelope smoothing: `lfp.power_smoothing_ms` (0 = off)
- Normalization: `lfp.normalize ∈ {'zscore','minmax','none'}`
- LFP history: `lfp.history_window_s` and optional `lfp.history_basis`
- Link: `lfp.link ∈ {'identity','log'}` — use `log` for strictly‑positive targets (e.g., power) if desired
- CV: `cv.k`, `cv.lambdas` — λ grid for smoothness strength
- Produced split: `produced_split_mode ∈ {'context','call_type'}`
- Streams: `perceived_addressed_window_s`, `perceived_overheard_silence_s`

Tips
- For high‑gamma power, keep `extract_band_power=true`, set the band to `[70,150]` Hz (or per‑array specifics), and consider a modest smoothing window (e.g., 25–50 ms) before resampling.
- Use a smaller `dt` for finer temporal alignment if signals and labels justify it (e.g., 10 ms).

---

## Shapes and Conventions

- `stim.t`: n_bins×1, centers, starting at 0, step `dt`
- LFP channel vector `lfp_ch`: n_bins×1 after `bin_lfp`
- `Xd.X`: n_bins×n_features (sparse), `Xd.y`: n_bins×1 (LFP), `Xd.colmap`: block metadata
- `D`: n_penalty_rows×n_features (sparse); `lambda` scales its effect
- Predictions are on the same grid as `stim.t`

---

## Practical Notes

- Out‑of‑range time bins (before/after the raw LFP window) are extrapolated to 0 during interpolation.
- History block excludes lag 0 by construction; set `lfp.history_window_s = [dt, L]`.
- To visualize LFP history weights today, either:
  - read them from `Xd.colmap.lfp_history.cols` in `wmap.w`, or
  - extend `src/model/unpack_params.m` and plots to include `lfp_history` like `spike_history`.
- Smoothness defaults penalize heard/produced kernel blocks; add `cfg.lfp_history=true` to also penalize the history block.

---

## What Gets Saved

`<outdir>/fit_results_lfp.mat` contains:
- `cfg`, `lfp_data` (meta), `events`
- `stim`, `lfp_binned`, `streams`, `states`
- Per‑channel `results(ch)` with `Xd`, `D`, `best_lambda`, `cvinfo`, `wmap`, `lfp_pred`, `lfp_actual`, `metrics`, `kernels`

Plots under `<outdir>/plots/`:
- `lfp_r2_by_channel.pdf`
- `lfp_traces_ch*.pdf`
- `lfp_kernels_summary_ch*.pdf` and `channel_*/kernels.pdf`

---

## Sanity Checks and Troubleshooting

- R² near 0 across channels: verify label alignment (`lfp_t0`), `dt`, and band selection; check that `produced_split_mode` matches your experiment.
- Unstable optimization: try `identity` link first; reduce λ grid if the minimum is at the edge; verify `D` is non‑empty for penalized blocks.
- No history in kernel plots: see the note above on `lfp_history` vs `spike_history` in the unpacker.

---

## References (Source Files)

- Runner: `scripts/run_fit_lfp_multichannel.m`
- I/O: `src/io/load_lfp.m`, `src/io/load_labels.m`
- Preprocess: `src/preprocess/build_timebase.m`, `src/preprocess/bin_lfp.m`, `src/preprocess/build_streams.m`, `src/preprocess/compute_states.m`, `src/preprocess/consolidate_twitter_bouts.m`
- Features: `src/features/assemble_design_matrix.m`, `src/features/smoothness_penalty.m`
- Model: `src/model/crossval_blocked.m`, `src/model/fit_glm_gauss.m`, `src/model/predict_lfp.m`, `src/model/unpack_params.m`
- Eval/plots: `src/eval/metrics_lfp.m`, `src/viz/plot_lfp_r2_by_channel.m`, `src/viz/plot_lfp_traces.m`, `src/viz/plot_lfp_kernels_summary.m`, `src/plot/plot_kernels.m`
