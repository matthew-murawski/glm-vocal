
# GLM for Natural Marmoset Conversation — Implementation Blueprint

**Document purpose.** This blueprint turns the SPEC into a concrete build plan with: (1) phased milestones, (2) progressively smaller task breakdowns, and (3) a series of **code-generation LLM prompts** to implement each step **test‑first**, with no orphaned code.

**Context.** We are re‑implementing a Li–Aoi–Miller–style GLM for Zhao & Wang (2023) natural conversational data. Use the repository layout and contracts defined in `SPEC.md` (provided). This blueprint assumes MATLAB (R2021b+) on macOS.

**Style constraints.** For MATLAB code, follow user’s comment style rules:
- section-by-section comments: lowercase, clear intention, not verbose.
- inline comments: minimal, only when necessary.
- use uppercase only for abbreviations (e.g., “ACC”).

---

## Guiding principles

- **Small, integrated steps.** Each step ends wired into a runnable script or test.
- **Design for sparsity.** Feature blocks & penalties must be sparse to scale.
- **Deterministic & validated.** Strong input validation, reproducible RNG, blocked CV.
- **Observability.** Plots + logs at key points (design matrix, CV curve, kernels).

---

## Phase overview (milestones)

**P0. Repo & scaffolding**  
Scaffold repo, config, and utilities; add CI-ready test harness (local).

**P1. I/O + validation**  
Load spikes & labels; normalize; validate shapes & ranges.

**P2. Time base & binning**  
Build uniform time grid; rasterize spikes to counts; add masks.

**P3. Streams & states**  
Rasterize produced/heard streams; compute conversational vs spontaneous state flags.

**P4. Features → design matrix**  
Build lagged kernel blocks (causal/symmetric); spike-history; assemble sparse `X`, `y`, col map.

**P5. Smoothness**  
Second-difference penalty blocks; per‑block λ support.

**P6. Fitting & CV**  
Poisson NLL + gradient; blocked K‑fold CV over λ; final MAP fit; prediction.

**P7. Evaluation & plots**  
Metrics, QC plots (kernels, rate vs spikes, CV curve, X preview).

**P8. End‑to‑end script**  
`run_fit_single_neuron.m` wires everything; save artifacts; logs.

**P9. Hardening**  
Error handling, collinearity checks, performance tuning, synthetic e2e test, mini‑real e2e.

**P10. Extensions (optional)**  
Per‑type kernels, nuisance streams, syllables, identity‑aware heard streams.

---

## Milestones → Iterative chunks

Each milestone is broken into **iterative chunks**; each chunk ends with a runnable test or demo.

### P0 — Repo & scaffolding
- **C0.1** Create folders & empty files per SPEC; add `README.md`, `docs/SPEC.md` link.
- **C0.2** Add a simple local test runner script (MATLAB) and `tests/` skeleton.
- **C0.3** Add `config/defaults.json` & `src/utils/struct_tools.m` helpers.

### P1 — I/O + validation
- **C1.1** `load_spikes.m` with unit tests on mock MAT files.
- **C1.2** `load_labels.m` reading MAT and Audacity txt; normalize to struct array.
- **C1.3** `validate_inputs.m` enforcing monotonic spikes, event interval sanity.

### P2 — Time base & binning
- **C2.1** `build_timebase.m` (t, dt, mask.good).
- **C2.2** `bin_spikes.m` exact histogram counts; edge-bin handling; tests.

### P3 — Streams & states
- **C3.1** `build_streams.m` constructing `heard_any` & `produced_any` from intervals; overlap-safe.
- **C3.2** `compute_states.m` rule‑based conversational vs spontaneous using response window; interval algebra & rasterization; tests.

### P4 — Features → design matrix
- **C4.1** `build_kernel_block.m` for causal streams.
- **C4.2** Extend to symmetric streams (produced).
- **C4.3** `build_history_block.m` from `sps` (exclude lag 0).
- **C4.4** `assemble_design_matrix.m` concatenating blocks + intercept + states → sparse `X`, `y`, `colmap`.
- **C4.5** `plot_design_matrix.m` quick visual QC.

### P5 — Smoothness
- **C5.1** `smoothness_penalty.m` second‑diff (D) per kernel block; no penalty on states/intercept.
- **C5.2** Per‑block λ handling; tests for block indexing/map.

### P6 — Fitting & CV
- **C6.1** `neglogli_poiss.m` + gradient (finite‑diff check).
- **C6.2** `fit_glm_map.m` (fminunc or LBFGS) with operator form for D' D; clipping safeguards.
- **C6.3** `crossval_blocked.m` K folds, λ grid, held‑out NLL per bin; choose λ; plot.
- **C6.4** `predict_rate.m`, `unpack_params.m`.

### P7 — Evaluation & plots
- **C7.1** `metrics.m` (NLL, pseudo‑R2, deviance explained).
- **C7.2** `plot_kernels.m`, `plot_rate_vs_spikes.m`, `plot_cv_curve.m`.
- **C7.3** `qc_session_summary.m` (one‑stop recap).

### P8 — End‑to‑end
- **C8.1** `scripts/run_fit_single_neuron.m` wiring + save artifacts.
- **C8.2** `scripts/demo_synthetic.m` (uses P9 synthetic generator).

### P9 — Hardening
- **C9.1** `check_collinearity.m` & warnings.
- **C9.2** Synthetic generator + e2e test (known kernels → Poisson spikes).
- **C9.3** Mini real‑data slice e2e test & smoke plots.
- **C9.4** Logging, reproducibility, config snapshot saved with results.

---

## Second breakdown — granular tasks (ready for TDD)

Below, each chunk is split into **small steps** that are independently testable.

### C0.1 Repo skeleton
1. Create directory tree (`src`, `tests`, `scripts`, `config`, `docs`, `plots`, `results`).
2. Add empty stubs matching SPEC signatures.
3. Add `README.md` quickstart.

### C0.2 Local test runner
1. Minimal `scripts/run_tests.m` that discovers and runs tests in `tests/`.
2. Ensure non‑zero exit on failure; print summary.

### C0.3 Config & utils
1. `config/defaults.json` with MVP values from SPEC.
2. `src/utils/struct_tools.m` (merge configs; get/set nested).

### C1.1 load_spikes
1. Implement MAT loader with fields validation.
2. Unit tests for well‑formed & malformed inputs.

### C1.2 load_labels
1. Implement MAT loader (`events` struct array).
2. Audacity `.txt` parser (on/off/label).
3. Normalizer: ensure `kind` ∈ {produced, perceived}, t_off≥t_on.
4. Unit tests for MAT & TXT paths; malformed rows.

### C1.3 validate_inputs
1. Spike monotonicity, finite, session bounds.
2. Label intervals in range; no NaNs.
3. Throw with identifiers; unit tests.

### C2.1 build_timebase
1. `dt`, T from spikes/labels; build `t`, `mask.good=true`.
2. Edge rounding; unit tests for length & last bin.

### C2.2 bin_spikes
1. `histcounts` vectorized; exact sum match to spikes.
2. Tests: empty spikes, spikes at edges, dense spikes.

### C3.1 build_streams
1. Interval→raster for produced & perceived; inclusive start, exclusive end.
2. Overlap tests; boundary tests.

### C3.2 compute_states
1. Response‑window rule; interval algebra (open/close windows).
2. Ensure mutual exclusivity; unit tests with synthetic dialogues.

### C4.1 build_kernel_block (causal)
1. Shifts via index arithmetic; sparse triplets.
2. Tests: correct lag alignment; off‑by‑one checks.

### C4.2 symmetric extension (produced)
1. Pre and post lags; center at event onset bin.
2. Tests: symmetry of lag indexing relative to onset.

### C4.3 build_history_block
1. Exclude lag 0; causal only; uses `sps`.
2. Tests with synthetic spike trains.

### C4.4 assemble_design_matrix
1. Concatenate blocks: intercept | heard | produced | states | history.
2. Build `colmap`; drop `~mask.good` rows.
3. Tests for column counts & map ranges.

### C4.5 plot_design_matrix
1. Non‑blocking preview (imshow of a column subset).
2. Smoke test only.

### C5.1 smoothness_penalty
1. Second‑difference matrix per kernel block.
2. Block‑diag assembly; no penalty on states/intercept.
3. Tests against expected finite difference stencils.

### C5.2 per‑block λ
1. Accept struct of λs; scale D rows accordingly.
2. Tests on index routing and scales.

### C6.1 neglogli_poiss
1. Return nll & grad; guard for exp overflow via clipping Xw.
2. Gradient check via finite differences.

### C6.2 fit_glm_map
1. fminunc wrapper; operator for D' D; early stopping.
2. Tests on tiny synthetic with known optimum.

### C6.3 crossval_blocked
1. Fold splitter; λ grid loop; aggregate NLL.
2. Tests: fold coverage, λ ordering, winner selection.

### C6.4 predict_rate & unpack_params
1. `mu=exp(X*w)`; unpack blocks via `colmap` into time vectors.
2. Unit tests for shapes & indices.

### C7.1 metrics
1. NLL, pseudo‑R2, deviance explained.
2. Tests on synthetic with known baselines.

### C7.2 plots
1. Kernels, rate vs spikes, CV curve; save to files.
2. Smoke tests (files exist).

### C7.3 qc_session_summary
1. Compose one‑stop figure & text summary; save.
2. Smoke test; basic content checks.

### C8.1 run_fit_single_neuron
1. Load config/spikes/labels; build stim/streams/states; assemble X; build D; CV; fit; predict; save artifacts; plots.
2. Integration test with synthetic/mini real.

### C8.2 demo_synthetic
1. Generate synthetic streams; known kernels; Poisson spikes; run full pipeline; check recovery metrics.

### C9.x hardening
1. Collinearity checker; warnings.
2. Logging; config snapshot; RNG control.
3. E2E regression test (fixed seed).

---

## Code‑generation LLM prompts (test‑driven).

> Use these sequentially. Each prompt is **self‑contained** and assumes prior steps exist. Keep MATLAB comments **lowercase**; use concise inline comments only when necessary.

### Prompt 0.1 — repo skeleton
```text
You are generating MATLAB files for a project named “glm-vocal”. Create the following folders and empty stubs matching these signatures. Do not write any business logic yet, only the function headers and minimal returns so tests can import them.

Folders:
- src/io, src/preprocess, src/features, src/model, src/eval, src/plot, src/utils
- tests/unit, tests/e2e, scripts, config, docs

Create stubs:
- src/io/load_spikes.m
- src/io/load_labels.m
- src/io/validate_inputs.m
- src/preprocess/build_timebase.m
- src/preprocess/bin_spikes.m
- src/preprocess/build_streams.m
- src/preprocess/compute_states.m
- src/features/build_kernel_block.m
- src/features/build_history_block.m
- src/features/assemble_design_matrix.m
- src/features/smoothness_penalty.m
- src/model/neglogli_poiss.m
- src/model/fit_glm_map.m
- src/model/crossval_blocked.m
- src/model/unpack_params.m
- src/model/predict_rate.m
- src/eval/metrics.m
- src/plot/plot_kernels.m
- src/plot/plot_rate_vs_spikes.m
- src/plot/plot_cv_curve.m
- src/plot/plot_design_matrix.m
- src/utils/struct_tools.m
- scripts/run_fit_single_neuron.m
- scripts/run_tests.m

Ensure each function compiles. Add minimal “help” lines describing purpose.
```

### Prompt 0.2 — config + utils
```text
Implement:
- config/defaults.json with MVP values:
  dt=0.01; heard_window_s=[0,0.4]; produced_window_s=[-0.5,0.5];
  history_window_s=[0.01,0.2]; cv.k=5; cv.lambdas=[0.1,0.3,1,3,10,30,100];
  state.response_window_s=5.0; state.state_window_s=5.0;
  optimizer.tol_fun=1e-6; optimizer.tol_grad=1e-5; optimizer.max_iter=200; seed=42

Implement src/utils/struct_tools.m with:
- function out = struct_merge(base, override)  % deep merge structs
- function v = struct_get(s, path, default)   % dot path get
- function s = struct_set(s, path, value)     % dot path set

Write unit tests in tests/unit/test_struct_tools.m.
```

### Prompt 1.1 — load_spikes (TDD)
```text
Implement src/io/load_spikes.m:
- Input: path to MAT with fields spike_times (double col vec), neuron_id, session_id
- Output: struct sp with same fields; validate types; spike_times sorted ascending

Unit tests in tests/unit/test_load_spikes.m:
- happy path MAT (create temp MAT in test)
- unsorted spikes -> function sorts and warns
- non-finite values -> error('glm:InvalidSpikeTimes', ...)
- missing fields -> error('glm:InvalidSpikesStruct', ...)
```

### Prompt 1.2 — load_labels (TDD)
```text
Implement src/io/load_labels.m:
- Accept a MAT file with events struct array OR an Audacity .txt
- Normalize to struct array with fields: kind ('produced'|'perceived'), t_on, t_off, label (string)
- Validate t_off>=t_on; kinds valid; finite numbers

Tests in tests/unit/test_load_labels.m cover MAT and TXT, malformed rows, and label normalization.
```

### Prompt 1.3 — validate_inputs (TDD)
```text
Implement src/io/validate_inputs.m that checks:
- spike_times finite & non-decreasing
- events finite, t_off>=t_on
- times within [0, T+buffer]

Throw using identifiers:
- glm:InvalidSpikeTimes
- glm:InvalidLabels

Write tests tests/unit/test_validate_inputs.m.
```

### Prompt 2.1 — build_timebase (TDD)
```text
Implement src/preprocess/build_timebase.m:
- Inputs: ev struct, sp struct (for T inference), dt
- Output stim struct with fields:
  stim.t (0:dt:T)', stim.dt=dt, stim.mask.good = true(size(t))

Tests test_build_timebase.m:
- correct length with exact last bin covering the last event/spike
- dt rounding behavior
```

### Prompt 2.2 — bin_spikes (TDD)
```text
Implement src/preprocess/bin_spikes.m:
- Inputs: spike_times, stim
- Output: sps, counts per bin via histcounts with edges [t; t(end)+dt]

Tests test_bin_spikes.m:
- sum(sps) equals numel(spike_times)
- spikes at bin edges
- empty spikes
```

### Prompt 3.1 — build_streams (TDD)
```text
Implement src/preprocess/build_streams.m:
- Inputs: ev, stim
- Outputs: streams struct with logical vectors heard_any, produced_any aligned to stim.t
- Definition: mark bin as 1 if interval [t_on,t_off) overlaps the bin

Tests test_build_streams.m:
- single interval
- overlapping intervals
- edge inclusion/exclusion conventions
```

### Prompt 3.2 — compute_states (TDD)
```text
Implement src/preprocess/compute_states.m:
- Inputs: ev, stim, stateCfg with response_window_s and state_window_s
- Output: states struct with logical vectors convo, spon (mutually exclusive)

Rule:
- on every perceived onset, open a response window of W_resp seconds
- if a produced onset occurs inside that window, mark conversational from perceived onset to produced offset
- otherwise mark perceived onset to perceived onset + state_window_s as conversational
- all other bins are spontaneous

Tests test_compute_states.m with synthetic timelines including overlaps and multiple perceived events.
```

### Prompt 4.1 — build_kernel_block (causal) (TDD)
```text
Implement src/features/build_kernel_block.m:
- Inputs: z (logical/double vector), stim, window_s=[0, L], mode='causal'
- Output: Xblk sparse [nT x Lbins], info with lag times

Test test_build_kernel_block.m:
- verify column k equals z shifted by k bins
- off-by-one checks, edges zero-padded
```

### Prompt 4.2 — build_kernel_block (symmetric) (TDD)
```text
Extend build_kernel_block to support mode='symmetric' with window_s=[-Lpre, Lpost].
- Center column corresponds to lag 0 at event onset
- Pre lags are negative, post lags positive

Update tests for symmetric alignment vs event onset indices.
```

### Prompt 4.3 — build_history_block (TDD)
```text
Implement src/features/build_history_block.m:
- Inputs: sps, stim, window_s=[dt, L_hist]
- Output: Xblk causal history without lag 0

Tests test_build_history_block.m on synthetic spike trains.
```

### Prompt 4.4 — assemble_design_matrix (TDD)
```text
Implement src/features/assemble_design_matrix.m:
- Inputs: streams, states, sps, cfg (windows), stim
- Build blocks:
  intercept, heard_any (causal), produced_any (symmetric), states (two scalar cols), spike_history (causal)
- Output Xd with fields X (sparse), y (sps), colmap

Tests test_assemble_design_matrix.m:
- correct column counts, order, and y alignment
- respects stim.mask.good if present
```

### Prompt 4.5 — plot_design_matrix
```text
Implement src/plot/plot_design_matrix.m for quick visualization of a subset of columns; save to file; smoke test.
```

### Prompt 5.1 — smoothness_penalty (TDD)
```text
Implement src/features/smoothness_penalty.m:
- Inputs: colmap, cfg (which blocks to penalize)
- Build second-difference matrices per kernel block; concatenate block-diagonal D; return D and Dmap

Tests test_smoothness_penalty.m validating stencils and index ranges; no penalty applied to states/intercept.
```

### Prompt 5.2 — per-block lambdas
```text
Allow lambda to be scalar or struct: lambda.heard, lambda.produced, lambda.history.
Implement scaling of D rows accordingly. Add tests for scaling behavior.
```

### Prompt 6.1 — neglogli_poiss (TDD)
```text
Implement src/model/neglogli_poiss.m:
- Inputs: w, X, y, and function handle L2op(w) that returns penalty value and gradient addition
- Output: nll, grad; clip Xw to [-50,50] to avoid overflow

Tests test_neglogli_poiss.m with finite-difference gradient check.
```

### Prompt 6.2 — fit_glm_map (TDD)
```text
Implement src/model/fit_glm_map.m:
- Use fminunc; objective = nll(X,w,y) + lambda*||D w||^2 via operator L2op
- Early stopping by tol_fun / tol_grad; return wmap.w, iterations, final objective

Tests test_fit_glm_map.m on a small synthetic X,y with known optimum; ensure convergence and positive rate predictions.
```

### Prompt 6.3 — crossval_blocked (TDD)
```text
Implement src/model/crossval_blocked.m:
- Inputs: fitfun handle, Xd, D, cvCfg with k and lambdas
- Split contiguous folds; for each lambda, fit on K-1 folds, compute test NLL per bin on held-out; pick best

Tests test_crossval_blocked.m verifying fold coverage, lambda selection on synthetic data.
```

### Prompt 6.4 — predict_rate & unpack_params (TDD)
```text
Implement src/model/predict_rate.m (mu=exp(X*w)) and src/model/unpack_params.m (slice w into named kernels using colmap).

Tests test_predict_and_unpack.m for shapes and map accuracy.
```

### Prompt 7.1 — metrics (TDD)
```text
Implement src/eval/metrics.m:
- Compute NLL, pseudo-R2, and deviance explained vs a constant-rate baseline

Tests test_metrics.m with known cases.
```

### Prompt 7.2 — plotting (smoke)
```text
Implement plot_kernels.m, plot_rate_vs_spikes.m, plot_cv_curve.m; each saves figures; smoke tests assert files exist.
```

### Prompt 7.3 — qc_session_summary
```text
Implement src/eval/qc_session_summary.m:
- Create and save a one-stop figure and a text/JSON summary with key stats and paths.
- Smoke test only.
```

### Prompt 8.1 — run_fit_single_neuron (integration)
```text
Implement scripts/run_fit_single_neuron.m:
- Load cfg, spikes, labels
- stim <- build_timebase; sps <- bin_spikes
- streams <- build_streams; states <- compute_states
- Xd <- assemble_design_matrix; [D, Dmap] <- smoothness_penalty
- [best_lambda, cvinfo] <- crossval_blocked
- [wmap, fitinfo] <- fit_glm_map
- mu <- predict_rate; kernels <- unpack_params
- metrics; plots; qc summary; save artifacts

Add tests/e2e/test_run_fit_single_neuron.m using synthetic data.
```

### Prompt 8.2 — demo_synthetic (integration)
```text
Implement scripts/demo_synthetic.m to generate synthetic streams, known kernels, and Poisson spikes, then call run_fit_single_neuron.m; assert recovery quality (pseudo-R2 > 0.05).
```

### Prompt 9.x — hardening
```text
Implement check_collinearity.m; integrate warnings into assemble_design_matrix or qc. Add logging and config snapshot saving. Create a regression e2e test with fixed seed and thresholded kernel L2 distance.
```

---

## Acceptance checklist (for each step)

- Unit/e2e tests exist and pass locally.
- New code is reachable from either a test or `run_fit_single_neuron.m`.
- No unreferenced files or dead code paths.
- Plots/outputs write into `results/` or `plots/`, not in source directories.

---

## Suggested implementation order (daily plan)

1. P0–P1 (Day 1): Skeleton, I/O, validation.
2. P2–P3 (Day 2): Time base, binning, streams, states.
3. P4 (Day 3): Feature blocks and design matrix.
4. P5–P6 (Day 4–5): Smoothness, NLL/grad, CV, fit.
5. P7–P8 (Day 6): Evaluation, plots, E2E script.
6. P9 (Day 7+): Hardening, synthetic & real-slice runs.
