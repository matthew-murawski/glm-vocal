
# GLM for Frontal Cortex in Natural Marmoset Conversation — Developer Specification

**Status:** v0.1 (MVP spec)  
**Owner:** you (Zhao Lab)  
**Language:** MATLAB
**Goal:** Rebuild a Li–Aoi–Miller–style GLM from scratch for Zhao & Wang (2023) naturalistic data (freely moving, multi-animal), keeping the core architecture but adapting event engineering to multi-speaker conversational scenes.

---

## 0) TL;DR (for the implementer)
- **Inputs:** per-neuron spike times; labeled produced calls (on/off, type optional); labeled perceived calls (on/off, type optional); optional noise/background regressors.  
- **Transform:** bin to fixed `dt`; make **regressor streams**: `heard_addressed` and `heard_overheard` (causal kernels split by conversational proximity), `produced_any` (symmetric kernel), `state_convo` / `state_spon` (scalars), plus **spike-history** (causal kernel). Build a **sparse design matrix** with lagged copies per stream.  
- **Model:** Poisson GLM with exponential link; **smoothness L2** per kernel; **blocked CV** to choose λ; fit MAP; unpack kernels; plot & evaluate.  
- **Outputs:** learned kernels (heard, produced, history), state coefficients, predicted rate, fit metrics, QC plots.  
- **MVP differences from Li:** “speaker” → “other-marmosets heard”; define conversational/spontaneous states by response-window rules; handle overlaps; start pooled across call types.

---

## 1) Scope & Requirements

### 1.1 MVP
- Single-session, single-neuron GLM fit.
- Regressors:
- **heard_addressed / heard_overheard:** binary time series marking perceived calls classified by subject-produced-call proximity (addressed: onset within 4 s of any subject production; overheard: at least 5 s of subject silence before and after). Both use **causal** kernels over `[0, L_heard]` ms.
  - **produced_any:** binary time series (1 when subject produces a call). **Symmetric** kernel over `[-L_prod_pre, +L_prod_post]` ms.
  - **state_convo, state_spon:** binary “context” flags with **scalar** coefficients (no temporal lags).
  - **spike_history:** causal kernel over `[dt, L_hist]` ms from the neuron’s own spike train.
- Model: Poisson GLM, log link (`exp` rate), L2 smoothness penalty on kernel blocks; no penalty on intercept or state scalars.
- Cross-validation: K-fold **blocked** (contiguous time folds); choose λ by held-out negative log-likelihood (NLL) per bin.
- Outputs: parameter vector `wmap`, unpacked kernels, predicted rate per bin, metrics, plots.
- Robust to overlaps (heard & produced may be 1 simultaneously).

### 1.2 Non-goals in MVP (roadmap)
- No per-identity regressors for specific marmosets (future).
- No per-call-type kernels initially (future split when counts suffice).
- No hierarchical/mixed-effects across neurons (future).
- No GLMnet/coordinate descent dependency; keep native MATLAB + sparse ops.

### 1.3 Data Assumptions
- Spike times in seconds, per neuron, relative to session start.
- Labeled events (produced & perceived) are (start_time, end_time, label_string). Labels optional for MVP. Times in seconds.
- Session duration T ≤ several hours; recording gaps permitted (we treat time uniformly; if gaps exist, mark as “bad” mask).

---

## 2) Repository Structure

```
glm-vocal/
├─ README.md
├─ docs/
│  └─ SPEC.md                      # this document
├─ config/
│  └─ defaults.json                # default hyperparameters (dt, kernel lengths, CV, etc.)
├─ src/
│  ├─ io/
│  │  ├─ load_spikes.m
│  │  ├─ load_labels.m             # reads produced/perceived labels into a normalized struct
│  │  └─ validate_inputs.m
│  ├─ preprocess/
│  │  ├─ build_timebase.m          # builds time grid, masks
│  │  ├─ bin_spikes.m
│  │  ├─ build_streams.m           # builds heard_addressed/overheard, produced_any, optional type-specific
│  │  └─ compute_states.m          # derives conversational/spontaneous flags from rules
│  ├─ features/
│  │  ├─ build_kernel_block.m      # makes lagged-columns for a binary stream
│  │  ├─ build_history_block.m
│  │  ├─ assemble_design_matrix.m  # concatenates blocks + intercept + states; returns sparse X
│  │  └─ smoothness_penalty.m      # builds block-diagonal second-diff matrices (D) and index map
│  ├─ model/
│  │  ├─ neglogli_poiss.m
│  │  ├─ fit_glm_map.m             # MAP with L2: min NLL + λ||Dw||^2
│  │  ├─ crossval_blocked.m
│  │  ├─ unpack_params.m           # slices wmap into named sub-vectors
│  │  └─ predict_rate.m
│  ├─ eval/
│  │  ├─ metrics.m                 # NLL, pseudo-R2, deviance explained
│  │  ├─ qc_session_summary.m
│  │  └─ check_collinearity.m
│  ├─ plot/
│  │  ├─ plot_kernels.m
│  │  ├─ plot_rate_vs_spikes.m
│  │  ├─ plot_cv_curve.m
│  │  └─ plot_design_matrix.m
│  └─ utils/
│     ├─ to_sparse_toeplitz.m
│     ├─ shifts_to_columns.m
│     ├─ mask_intervals.m
│     └─ struct_tools.m
├─ tests/
│  ├─ unit/
│  └─ e2e/
└─ scripts/
   ├─ run_fit_single_neuron.m      # CLI/script entry point for a session
   └─ demo_synthetic.m
```

---

## 3) Data Handling & Formats

### 3.1 Spike Times
```matlab
% spikes.mat
spike_times: double column vector [n_spikes x 1], seconds from session start
neuron_id:   char or string scalar
session_id:  char or string scalar
```

### 3.2 Labels (Produced, Perceived)
Accept one or multiple label files; combine into a normalized table/struct.

```matlab
% labels.mat (or parsed from Audacity txt)
events: struct array with fields:
  .kind        % 'produced' | 'perceived'
  .t_on        % double, seconds
  .t_off       % double, seconds (>= t_on)
  .label       % call type or free text (optional for MVP)
  .quality     % optional flag, e.g., 'ok'|'uncertain'|'noise'
session_id: string
subject_id: string
```

### 3.3 Time Base
- Choose `dt` (default: 0.01 s).
- `t = (0:dt:T)'` where `T` is the max of spikes/labels.  
- Maintain a `mask.good` vector (logical) to drop bad/unknown periods if needed.

```matlab
stim.t: [nT x 1] double
stim.dt: scalar double
stim.mask.good: [nT x 1] logical (default all true)
```

### 3.4 Regressor Streams
From `events` to binary time series aligned to `stim.t`:
- `streams.heard_addressed` / `streams.heard_overheard` (1 at perceived-call onset bins classified as addressed vs overheard)
- `streams.produced_any` (1 if any subject call overlaps the bin interval)
- Optional: per-type versions later (e.g., `heard_phee`, `produced_twitter`), **not in MVP**.

### 3.5 State Flags (Conversational vs Spontaneous)
**Rule (MVP):**
- When a **perceived** call occurs, open a **response window** of `W_resp` seconds (default 5 s).
- If a **produced** call occurs with onset inside that window, mark a **conversational state** from the perceived onset through the produced offset (or through a fixed `W_state` secs after perceived onset if no produced call occurs).  
- Bins not in conversational are **spontaneous**.
- Export `states.convo` and `states.spon` (logical vectors). Ensure `states.convo & states.spon` is never true.

Implementation tip: build interval lists, do interval algebra, then rasterize to the time grid. Keep it deterministic and unit-tested.

### 3.6 Binning Spikes
- `sps = histcounts(spike_times, [t; t(end)+dt])'` to get counts per bin.
- Sanity check: `sum(sps) == numel(spike_times)`.

---

## 4) Feature Engineering → Design Matrix

### 4.1 Kernel Definitions (MVP defaults)
- **heard_addressed** / **heard_overheard** (causal): window `[0, 0.4]` s → `L_heard = 0.4/dt + 1` columns (including lag 0) for each stream.
- **produced_any** (symmetric): window `[-0.5, +0.5]` s → `L_pre = 0.5/dt`, `L_post = 0.5/dt` columns; center at t=0.
- **spike_history** (causal): window `[dt, 0.2]` s → `L_hist = 0.2/dt` columns (exclude lag 0).
- **states**: `state_convo`, `state_spon` → **two scalar columns** (no lags).
- **intercept**: single column of ones.

### 4.2 Column Construction
For each binary stream `z` (heard/prod):
- **causal (heard/history):** build columns by shifting `z` backward in time: `z(t - k*dt)`, `k=0..L-1` (with zeros before start).
- **symmetric (produced):** build pre and post columns: `z(t + k*dt)` for `k<0` (pre) and `z(t - k*dt)` for `k≥0` (post). Center column corresponds to t=0 (onset bin).

Use **sparse** assembly:
- Prefer preallocating `(nT, n_cols)` with `spalloc` and filling with index triplets.
- Utility `shifts_to_columns.m` builds the block given `(z, lags, causal|symmetric)`.

### 4.3 Design Matrix Layout (column order)
```
[ intercept | heard_addressed_block | heard_overheard_block | produced_any_block | state_convo | state_spon | spike_hist_block ]
```
Return:
```matlab
Xd.X         % sparse [nT x p]
Xd.colmap    % struct with fields, each = [start_idx end_idx]
Xd.y         % = sps
Xd.info      % dt, windows, etc.
```

### 4.4 Masks & Censoring
- If `mask.good` exists, **drop** those rows consistently from `X` and `y`.
- Optionally censor bins around labeling gaps with `mask_intervals.m`.

---

## 5) Model & Fitting

### 5.1 Poisson GLM
- Rate per bin: `mu = exp(X * w)`.
- Negative log-likelihood per bin: `nll = sum(mu - y .* log(mu) + log(y!))` (ignore constant term).

### 5.2 Smoothness Penalty
- Build **second-difference** penalty for each **kernel block** (heard, produced, history). No penalty on **intercept** or **state** columns.
- Let `D` be block-diagonal; objective: `J(w) = NLL(w) + λ * ||D w||^2`.
- Support **per-block λ** (e.g., stronger smoothing for produced vs history) via scaling rows of `D` or separate hyperparameters.

### 5.3 Optimization
- Use MATLAB’s `fminunc` or custom L-BFGS; provide gradient:
  - `grad = X'*(exp(Xw) - y) + 2*λ*(D'*D*w)`.
- Initialize `w = 0`.
- Convergence: relative change in `J` < 1e-6 or gradient norm < 1e-5, max 200 iters.
- **Numerics:** compute `Xw` efficiently; keep everything sparse; avoid forming `D'*D` explicitly (apply as operator).

### 5.4 Cross-Validation (Blocked)
- Split time into K contiguous folds (default K=5).
- For each λ in a grid (e.g., logspace(-2, 3, 8)):
  - Fit on K-1 folds; evaluate **test NLL per bin** on the held-out fold (no refit of λ per fold).
- Aggregate mean test NLL; choose λ with minimum (± 1-SE rule optional). Keep per-block λs if configured.

### 5.5 Final Fit
- Refit on full dataset using chosen λ(s).
- Return `wmap`, along with Hessian diag approximation (for SEs) if requested.

---

## 6) Outputs & Artifacts

- `results/wmap.mat`: struct with fields:
  ```matlab
  wmap.w              % [p x 1]
  wmap.colmap         % same col map as Xd
  wmap.dt             % bin size
  wmap.lambdas        % chosen λ (or per-block λs)
  wmap.metrics        % train/test NLL, pseudo-R2, etc.
  wmap.session_id, wmap.neuron_id
  ```
- `results/kernels.mat`: unpacked kernels as time vectors:
  - `heard_kernel.t, heard_kernel.w`
  - `produced_kernel.t, produced_kernel.w` (pre/post included)
  - `history_kernel.t, history_kernel.w`
  - `state_coeffs = [beta_convo, beta_spon]`
- `plots/`:
  - `cv_curve.pdf`
  - `kernels.pdf` (multi-panel)
  - `rate_vs_spikes.pdf`
  - `design_matrix.pdf`

---

## 7) Configuration & Hyperparameters

`config/defaults.json` (example):
```json
{
  "dt": 0.01,
  "heard_window_s": [0.0, 0.4],
  "produced_window_s": [-0.5, 0.5],
  "history_window_s": [0.01, 0.2],
  "cv": { "k": 5, "lambdas": [0.1, 0.3, 1, 3, 10, 30, 100] },
  "state": { "response_window_s": 5.0, "state_window_s": 5.0 },
  "optimizer": { "tol_fun": 1e-6, "tol_grad": 1e-5, "max_iter": 200 },
  "seed": 42
}
```
A MATLAB helper `load_config.m` can read JSON and merge with overrides.

---

## 8) Error Handling & Validation

- **Input validation (`validate_inputs.m`):**
  - spikes must be increasing, finite; warn on duplicates.
  - labels must have `t_off >= t_on` and lie within `[0, T+buffer]`.
  - throw with identifier strings, e.g., `glm:InvalidSpikeTimes`, `glm:InvalidLabels`.
- **Empty/rare-event protection:**
  - If a stream is all zeros → drop that block and log a warning.
  - If event count < min_events (default 30) for a block → warn and consider pooling.
- **Collinearity diagnostics:**
  - Compute `diag(X'X)` norms and pairwise correlations among column blocks; warn if |ρ| > 0.95 between blocks.
- **Numerical stability:**
  - Clip `Xw` to [-50, 50] for `exp` to avoid overflow; track if clipping occurred.
- **Masking consistency:**
  - Ensure `X` and `y` undergo identical row drops; assert sizes match.
- **CV splits:**
  - Require each fold to contain at least N events per kernel block; else adjust K or merge folds.

---

## 9) Testing Plan

### 9.1 Unit Tests (`tests/unit`)
- `test_build_timebase.m`: correct length, dt, masks.
- `test_bin_spikes.m`: exact histogram match; edge bins.
- `test_build_streams.m`: interval → raster logic; overlaps; boundaries.
- `test_compute_states.m`: response-window rules on synthetic timelines.
- `test_build_kernel_block.m`: lag indexing; causal vs symmetric alignment.
- `test_smoothness_penalty.m`: second-diff stencil; block-diagonal assembly; no penalty on scalar columns.
- `test_neglogli_poiss.m`: values vs finite-difference gradient checks.
- `test_fit_glm_map.m`: convergence on small synthetic; gradient correctness.
- `test_crossval_blocked.m`: fold boundaries; λ selection monotonicity on a synthetic with known optimum.
- `test_unpack_params.m`: column-map slicing correct.
- `test_predict_rate.m`: shapes and positivity; invariance to intercept shifts.

### 9.2 End-to-End (`tests/e2e`)
- **Synthetic conversation**: generate perceived → produced with known GLM weights; simulate spikes via Poisson; verify recovery (up to smoothing) of kernel shapes and positive pseudo-R2.
- **Mini real session**: run on a 5–10 minute slice; ensure no errors; outputs created; CV curve non-degenerate.

### 9.3 Regression Tests
- Pin a small synthetic seed and check kernel L2 distance < ε after code changes.

---

## 10) Developer Ergonomics

- **Logging:** simple logger printing major steps (start/end times, n events, λ grid, best λ, NLLs).
- **Reproducibility:** set global RNG; save config snapshot with results; record git commit hash if available.
- **Speed & memory:** sparse matrices everywhere; build design matrix in chunks if needed; avoid dense `D'*D` via function handle application.
- **Style:** follow MATLAB style; for in-file code comments conform to user’s preference (lowercase section comments; concise inline only when needed).

---

## 11) Example Usage

```matlab
% scripts/run_fit_single_neuron.m
cfg = load_config("config/defaults.json");
sp = load_spikes("data/sessionA/neuron001_spikes.mat");
ev = load_labels("data/sessionA/labels.mat");

stim = build_timebase(ev, sp, cfg.dt);        % t, dt, mask
sps  = bin_spikes(sp.spike_times, stim);

streams = build_streams(ev, stim);            % heard_addressed/overheard, produced_any
states  = compute_states(ev, stim, cfg.state);

Xd = assemble_design_matrix(streams, states, sps, cfg);  % X, y, colmap

[D, Dmap] = smoothness_penalty(Xd.colmap, cfg);          % block-diag D

[best_lambda, cvinfo] = crossval_blocked(@fit_glm_map, Xd, D, cfg.cv);

[wmap, fitinfo] = fit_glm_map(Xd, D, best_lambda, cfg.optimizer);

rate = predict_rate(Xd.X, wmap.w);

kernels = unpack_params(wmap, Xd.colmap, cfg, stim);

qc_session_summary(Xd, wmap, rate, cvinfo, "results/");
```

---

## 12) Roadmap Extensions

- **Per-type kernels** (produced/heard) when counts permit.
- **Identity-aware heard regressors** once caller ID is annotated.
- **Nuisance covariates:** broadband level, noise events, slow drifts (splines).
- **Syllable-level regressors** for twitter calls.
- **Hierarchical GLM** across neurons or sessions; random effects on state coefficients.
- **Bayesian alternatives** with priors on smoothness; ARD over columns.
- **Model comparison**: AIC/BIC, held-out likelihood, nested tests for adding regressors.

---

## 13) Acceptance Criteria (MVP)

- Runs on one real session & neuron without errors.
- Produces: `wmap.mat`, `kernels.mat`, and plots in `results/` & `plots/`.
- Cross-validation completes; non-trivial λ selected; predicted rate > 0, finite.
- Kernels show sensible shapes (heard causal bump; produced symmetry around 0; decaying history).
- Pseudo-R2 on held-out > 0.05 (tunable threshold for MVP).
- Unit tests: ≥ 90% pass; E2E synthetic recovers ground-truth qualitatively.

---

## 14) File/Function Stubs (names + signatures)

```matlab
% src/io/load_spikes.m
function sp = load_spikes(path)

% src/io/load_labels.m
function ev = load_labels(pathOrPaths)

% src/io/validate_inputs.m
function validate_inputs(sp, ev)

% src/preprocess/build_timebase.m
function stim = build_timebase(ev, sp, dt)

% src/preprocess/bin_spikes.m
function sps = bin_spikes(spike_times, stim)

% src/preprocess/build_streams.m
function streams = build_streams(ev, stim)

% src/preprocess/compute_states.m
function states = compute_states(ev, stim, stateCfg)

% src/features/build_kernel_block.m
function [Xblk, info] = build_kernel_block(z, stim, window_s, mode)

% src/features/build_history_block.m
function [Xblk, info] = build_history_block(sps, stim, window_s)

% src/features/assemble_design_matrix.m
function Xd = assemble_design_matrix(streams, states, sps, cfg)

% src/features/smoothness_penalty.m
function [D, Dmap] = smoothness_penalty(colmap, cfg)

% src/model/neglogli_poiss.m
function [nll, grad] = neglogli_poiss(w, X, y)

% src/model/fit_glm_map.m
function [wmap, fitinfo] = fit_glm_map(Xd, D, lambda, optCfg)

% src/model/crossval_blocked.m
function [best_lambda, cvinfo] = crossval_blocked(fitfun, Xd, D, cvCfg)

% src/model/unpack_params.m
function kernels = unpack_params(wmap, colmap, cfg, stim)

% src/model/predict_rate.m
function mu = predict_rate(X, w)

% src/eval/metrics.m
function out = metrics(y, mu)

% src/plot/plot_kernels.m
function plot_kernels(kernels, outdir)

% src/plot/plot_rate_vs_spikes.m
function plot_rate_vs_spikes(stim, y, mu, outdir)

% src/plot/plot_cv_curve.m
function plot_cv_curve(cvinfo, outdir)

% src/plot/plot_design_matrix.m
function plot_design_matrix(Xd, outdir)
```

---

### glossary
- **lag:** how far a predictor from the past (or future, for symmetric) influences the current bin.
- **regressor stream:** a binary (or continuous) time series aligned to the bin grid that we expand into many lagged columns to learn a temporal **kernel**.
- **causal vs symmetric:** causal uses only present/past; symmetric adds pre- and post-event lags around an event time (e.g., produced onset).

---

**end of spec**
