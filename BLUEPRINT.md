# BLUEPRINT.md — glm-vocal (v1)
GLM for frontal-cortex spiking in natural vocal behavior, packaged as a reusable MATLAB library.

## Goal
- **Model**: Poisson GLM on **10 ms** bins (log link), predicting spike counts for a **single unit** over the **entire session**.
- **Covariates**:
  1) **Spike-history** (0–200 ms, 6 raised-cosine bases + 0–2 ms refractory boxcar).
  2) **Event kernels** for **produced** and **heard** call **onsets**, **interacting** with **conversational state** (conversational vs independent). Window −0.5…+1.5 s with 8 raised-cosine bases per kernel.
  3) **State term**: conversational-bout **boxcar** (slow covariate).
- **Conversational bouts**: runs where **inter-event gap < 4 s** across any calls (produced or heard).
- **Regularization**: **ridge (L2)** with **first-difference smoothing** within kernel/history groups.
- **Evaluation**: **5-fold blocked CV** for λ selection + final **20% hold-out**.

## Repository Layout
```
glm-vocal/
├─ src/+vglm/           % MATLAB package namespace (all library code lives here)
├─ tests/unit/          % unit tests (small, fast, synthetic)
├─ tests/e2e/           % end-to-end tests (short real/sim snippets)
├─ scripts/             % orchestration scripts (e.g., run_glm_single_unit.m)
├─ results/             % CSV/MAT artifacts
├─ figures/             % plots (filters, CV curves, diagnostics)
├─ docs/                % optional docs/notes
└─ README.md
```
**Comment style**: lower-case section comments; concise inline comments only when necessary.

## Data Flow (single session, single unit)
1) **Inputs**: spike times (s), produced/heard call labels (onset/offset/type/actor), session duration.
2) **Timebase**: build 10 ms grid → `t` with `vglm.build_timebase`.
3) **Spike binning**: `vglm.bin_spikes` → counts per bin `y`.
4) **Event cleanup**: `vglm.clean_events` (merge short gaps; drop ultra-short; dedupe).
5) **Bouts**: `vglm.compute_bouts` with gap < 4 s → conversational intervals; `vglm.bouts_to_state` → `state_conv(t)`.
6) **Event streams**: `vglm.events_to_streams` → four onset streams by state (prod/heard × conv/ind) + `state_conv`.
7) **Bases**: `vglm.make_basis_raised_cosine` + `vglm.sample_basis_on_grid` for history & event kernels; refractory boxcar for history.
8) **Design matrix**: `vglm.build_design_matrix` concatenates blocks (history, event kernels, state, intercept).
9) **Penalty**: `vglm.build_L2_matrix` builds block-wise first-difference matrices (intercept unpenalized).
10) **Fit**: `vglm.fit_poiss_glm_ridge` uses Poisson GLM with ridge; choose λ via blocked CV (`vglm.split_blocks`); refit; evaluate on final 20% hold-out.
11) **Artifacts**: metrics table, β coefficients, figures, and config MAT/CSV.

## Validations & Error Handling
- Assert bin width == 0.01 s and time coverage (start/end aligned).
- Labels: finite; onset < offset; monotonic; log drops/merges (counts and reasons).
- Overlap rules (same-actor high-overlap → merge) deterministic and documented.
- Design integrity: no NaN/Inf; blocks non-empty unless expected.
- Optimizer guardrails: retries from small-norm init if not converged; fail with clear report.
- CV integrity: contiguous folds; strict separation of final 20% hold-out.

## Metrics & Outputs
- **Selection**: mean test NLL vs λ (CV curve, chosen λ).
- **Goodness-of-fit**: hold-out NLL, pseudo-R², bits/spike vs history-only baseline.
- **Interpretability**: history and event kernels vs lag; state effect plot.
- **Reproducibility**: MAT (`beta`, bases, design meta), CSV summary, PNG plots in `figures/`.

---

# Iterative Build Plan → Chunks → Micro-Steps

## Phase A — Foundations (IO, timebase, binning)
**Chunk A0 — repo scaffold**
- Create folders: `src/+vglm/`, `tests/unit/`, `tests/e2e/`, `scripts/`, `results/`, `figures/`, `docs/` (optional).
- `scripts/startup.m`: add `src` to path; assert/create `results/`, `figures/`.
- `README.md`: overview + how to run tests:
  ```
  matlab -batch "addpath('scripts'); startup; results = runtests('tests','IncludeSubfolders',true); disp(results);"
  ```

**Chunk A1 — timebase & spike binning**
- A1.1 `vglm.build_timebase(session_dur, dt)`
- A1.2 `vglm.bin_spikes(spike_times, t)`
- A1.3 Unit tests (edges, empty, histogram parity)
- A1.4 Mini script to exercise A1

## Phase B — Event Preparation (cleanup & bouts)
**Chunk B1 — event cleanup**
- B1.1 `vglm.clean_events(tbl_in, params)` (merge gaps, min duration, dedupe/overlap)
- B1.2 Unit tests (merge, drop, overlap)
- B1.3 Mini script (before/after counts)

**Chunk B2 — conversational bouts**
- B2.1 `vglm.compute_bouts(prod_tbl, heard_tbl, gap_s)`
- B2.2 `vglm.bouts_to_state(bouts, t)`
- B2.3 Unit tests (toy timelines → bouts/state)

## Phase C — Bases & Design Blocks
**Chunk C1 — raised-cosine bases**
- C1.1 `vglm.make_basis_raised_cosine(twin, nb, densify_near0)`
- C1.2 `vglm.sample_basis_on_grid(spec, grid_s)`
- C1.3 Unit tests (shape, coverage, edges)

**Chunk C2 — event streams & kernel blocks**
- C2.1 `vglm.events_to_streams(prod_tbl, heard_tbl, bouts, t)` → 4 onset streams + `state_conv`
- C2.2 `vglm.convolve_events_with_basis(impulses, spec, t, event_win)`
- C2.3 Unit tests (single impulse → basis; superposition; state routing)

**Chunk C3 — spike-history block**
- C3.1 `vglm.build_history_block(y, dt, hist_win, nb, refractory)`
- C3.2 Unit tests (self-bin exclusion; refractory activity; shapes)

**Chunk C4 — design assembly**
- C4.1 `vglm.build_design_matrix(blocks)` → `X, meta` (group/name/lag, penalty_group_id)
- C4.2 Unit tests (counts, meta integrity, no NaN/Inf)

## Phase D — Penalty, Objective, Fitting, CV
**Chunk D1 — ridge penalty with smoothing**
- D1.1 `vglm.build_L2_matrix(meta)` → blockwise first-difference; intercept unpenalized
- D1.2 Unit tests (Δ structure; zero for intercept; dims)

**Chunk D2 — Poisson objective**
- D2.1 `vglm.poiss_glm_neglogpost(beta, X, y, dt, D, lambda)`
- D2.2 Finite-difference tests for gradient/Hessian

**Chunk D3 — fitter with blocked CV & final refit**
- D3.1 `vglm.split_blocks(T, kfold, test_frac)` → contiguous folds + final 20% hold-out
- D3.2 `vglm.fit_poiss_glm_ridge(X, y, dt, D, lambdas, cv)` → chosen λ, β̂, CV table, holdout metrics
- D3.3 Unit tests (synthetic truth; selection sanity)

## Phase E — Evaluation, Plots, Wiring
**Chunk E1 — evaluation & plots**
- E1.1 `vglm.evaluate_glm(X, y, dt, beta, baseline_rate)` → NLL, pseudo-R², bits/spike
- E1.2 `vglm.plot_filters(beta, meta, basis_info, out_path)`
- E1.3 `vglm.plot_cv_curves(cv_table, lambdas, out_path)`
- E1.4 Unit tests (smoke tests; PNGs saved)

**Chunk E2 — end-to-end script & e2e test**
- E2.1 `scripts/run_glm_single_unit.m` — wires everything for one session & unit
- E2.2 `tests/e2e/test_run_glm_single_unit.m` — 60–120 s snippet; asserts: finite λ; hold-out NLL < constant; sensible history
- E2.3 Save artifacts to `results/` and `figures/`

## Default Parameters (configurable)
```matlab
params = struct( ...
  'dt', 0.01, ...
  'hist_win', [0, 0.2], ...
  'n_hist', 6, ...
  'refractory', [0, 0.002], ...
  'event_win', [-0.5, 1.5], ...
  'n_event', 8, ...
  'bout_gap', 4.0, ...
  'merge_gap', 0.10, ...
  'min_dur', 0.05, ...
  'overlap_merge', 0.80, ...
  'lambdas', logspace(-2, 4, 20), ...
  'kfold', 5, ...
  'test_frac', 0.20 ...
);
```

---

# LLM Prompts (Test-Driven, No Orphans)

```text
[Prompt A0] — repo scaffold for glm-vocal
Goal: Initialize repo structure and MATLAB pathing.
Tasks:
- Create folders: src/+vglm/, tests/unit/, tests/e2e/, scripts/, results/, figures/, docs/.
- scripts/startup.m: add src to path; assert/create results/ and figures/.
- README.md: overview and how to run tests:
  matlab -batch "addpath('scripts'); startup; results = runtests('tests','IncludeSubfolders',true); disp(results);"
Style: user's comment rules.
```

```text
[Prompt A1] — timebase and spike binning
Goal: Foundational timebase and spike binning with tests.
Tasks:
1) src/+vglm/build_timebase.m — function t = vglm.build_timebase(session_dur, dt)
2) src/+vglm/bin_spikes.m — function y = vglm.bin_spikes(spike_times, t)
3) tests/unit/test_timebase_and_binning.m — edges, empty, histogram parity
4) scripts/dev_check_A1.m — exercise A1
```

```text
[Prompt B1] — event cleanup utilities
Goal: Normalize produced/heard label tables.
Tasks:
1) src/+vglm/clean_events.m — [tbl_out, loginfo] = vglm.clean_events(tbl_in, params)
2) tests/unit/test_clean_events.m — crafted merge/drop/overlap cases
3) scripts/dev_check_B1.m — before/after counts
```

```text
[Prompt B2] — conversational bouts and state boxcar
Goal: 4-s gap–based conversational bouts → state time series.
Tasks:
1) src/+vglm/compute_bouts.m — bouts = vglm.compute_bouts(prod_tbl, heard_tbl, gap_s)
2) src/+vglm/bouts_to_state.m — state = vglm.bouts_to_state(bouts, t)
3) tests/unit/test_bouts.m — toy timelines
```

```text
[Prompt C1] — raised-cosine bases
Goal: Basis generator and sampler.
Tasks:
1) src/+vglm/make_basis_raised_cosine.m — spec struct (centers, widths, twin, nb, densify_near0)
2) src/+vglm/sample_basis_on_grid.m — Phi = vglm.sample_basis_on_grid(spec, grid_s)
3) tests/unit/test_bases.m — shape, coverage, edges
```

```text
[Prompt C2] — event streams & kernel blocks
Goal: Four onset streams by state; kernel blocks via convolution.
Tasks:
1) src/+vglm/events_to_streams.m — streams + state_conv
2) src/+vglm/convolve_events_with_basis.m — Xblk from impulses × basis over event_win
3) tests/unit/test_event_streams_and_blocks.m — impulse→basis; superposition; state routing
```

```text
[Prompt C3] — spike-history block
Goal: 6 raised-cosines + refractory boxcar.
Tasks:
1) src/+vglm/build_history_block.m — Xhist = vglm.build_history_block(y, dt, hist_win, nb, refractory)
2) tests/unit/test_history_block.m — self-bin exclusion; refractory behavior
```

```text
[Prompt C4] — design matrix assembly & metadata
Goal: Concatenate blocks and tag columns.
Tasks:
1) src/+vglm/build_design_matrix.m — [X, meta] = vglm.build_design_matrix(blocks)
2) tests/unit/test_design_matrix.m — counts, meta integrity, no NaN/Inf
```

```text
[Prompt D1] — ridge penalty with smoothing
Goal: First-difference smoothing within groups; intercept unpenalized.
Tasks:
1) src/+vglm/build_L2_matrix.m — D = vglm.build_L2_matrix(meta)
2) tests/unit/test_penalty_matrix.m — Δ structure; intercept zero; dims
```

```text
[Prompt D2] — Poisson objective, gradient, Hessian
Goal: Stable NLL + ridge with derivatives.
Tasks:
1) src/+vglm/poiss_glm_neglogpost.m — [f,g,H] = vglm.poiss_glm_neglogpost(...)
2) tests/unit/test_neglogpost.m — finite-difference checks
```

```text
[Prompt D3] — fitter with blocked CV and final refit
Goal: Fit λ grid via contiguous CV; refit best; evaluate hold-out.
Tasks:
1) src/+vglm/split_blocks.m — cv = vglm.split_blocks(T, kfold, test_frac)
2) src/+vglm/fit_poiss_glm_ridge.m — model = vglm.fit_poiss_glm_ridge(...)
3) tests/unit/test_fit_cv.m — synthetic truth; selection sanity
```

```text
[Prompt E1] — evaluation & plotting
Goal: Compute metrics and produce figures.
Tasks:
1) src/+vglm/evaluate_glm.m — metrics = vglm.evaluate_glm(...)
2) src/+vglm/plot_filters.m — vglm.plot_filters(...)
3) src/+vglm/plot_cv_curves.m — vglm.plot_cv_curves(...)
4) tests/unit/test_eval_and_plots.m — smoke tests; PNGs saved
```

```text
[Prompt E2] — end-to-end wiring & e2e test (single unit)
Goal: Full pipeline from a script.
Tasks:
1) scripts/run_glm_single_unit.m — wires everything for one session & unit
2) tests/e2e/test_run_glm_single_unit.m — 60–120 s snippet; asserts on λ & NLL
3) Save CSV/MAT to results/; PNG to figures/; concise console summary
```

---

# Notes
- No orphaned utilities: each prompt yields runnable code and tests.
- Risky math (objective/grad/Hess) is isolated and verified before integration.
- CV and plotting come only after core correctness is established.
- Ready for extension to multi-unit models, offset kernels, and richer state variables later.
