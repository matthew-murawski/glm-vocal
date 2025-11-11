# High Gamma GLM — Implementation Checklist (Steps 1–14)

> Sectioned by **phase → step**. Each item is a concrete action or verification you can check off as you implement, test, and demo.

---

## Phase 0 — Foundations (Steps 1–2)

### Step 1 — Load and Validate High Gamma Data
**What to implement**
- [X] `src/io/load_high_gamma.m`
  - [X] Load MAT with `lfp_data.lfp_data [nSamples×nChannels]`, `lfp_data.lfp_fs`, `lfp_data.lfp_t0`, `lfp_data.channel_ids`, optional `metadata.timestamps`
  - [X] Support channel strategies: single `channel_id` **or** reduction across channels (`mean`, `median`)
  - [X] Return struct with fields: `power [nSamples×1]`, `fs`, `t [nSamples×1]`, `session_id` (if available), `channel_id` (selected or `'reduced'`)
  - [X] Time handling: prefer `metadata.timestamps`; else build `t = lfp_t0 + (0:n-1)'/lfp_fs`
  - [X] Ensure `power` is non-negative (clip/warn if negatives)
- [X] `src/io/validate_inputs_hg.m`
  - [X] Confirm `power` numeric vector, all finite
  - [X] Warn and set to zero if any negative values
  - [X] Error if all zeros **or** `std(power) < eps`
  - [X] Validate `fs > 0`
  - [X] Verify `t` strictly monotonically increasing
  - [X] Return validation struct (warnings/errors)
- [X] `tests/synthetic/generate_simple_hg_session.m`
  - [X] Create flat constant power trace with required fields populated
  - [X] Save as `tests/data/synthetic_simple.mat`

**Tests to write**
- [X] `tests/unit/test_load_high_gamma.m`
  - [X] Loads valid MAT successfully
  - [X] Single-channel selection works
  - [X] Multi-channel reduction (`mean`, `median`) works
  - [X] Time vector from `metadata.timestamps` branch
  - [X] Time vector from `lfp_t0`/`lfp_fs` branch
  - [X] Missing fields → informative errors
  - [X] Invalid types → informative errors
  - [X] Negative values trigger warning & clipping
  - [X] All-zero data → error
- [X] `tests/e2e/test_io_pipeline.m`
  - [X] Generate synthetic
  - [X] Load and validate
  - [X] Output struct fields and basic props verified

**Demo**
- [X] `scripts/dev/demo_step1_io.m`
  - [X] Verbose prints (dims, range, fs)
  - [X] Show head/tail samples
  - [X] Simple timeseries plot
  - [X] Print "✓ Step 1 Complete: I/O working"

**Success criteria**
- [X] All unit tests pass
- [X] Demo runs w/o errors
- [X] Single/reduced channel both supported
- [X] Correct time vector construction
- [X] Validation catches common errors

---

### Step 2 — Bin High Gamma to Time Base
**What to implement**
- [X] `src/preprocess/bin_high_gamma.m`
  - [X] Signature: `binned = bin_high_gamma(hg, stim, cfg)`
  - [X] Methods: `mean` | `median` | `rms` | `max` (from `cfg.preprocess.hg_binning_method`)
  - [X] Aggregate samples per `stim.t` bin
  - [X] Handle empty bins via linear interpolation (warn, count empties)
  - [X] If `hg.fs < 1/dt`, warn about upsampling; prefer `dt ≥ 1/fs`
  - [X] Ensure exact `length(stim.t)` and non-negative output
  - [X] Return fields: `power`, `t`, `method`, `n_empty_bins`, `warnings`
- [X] Helpers in same file
  - [X] `aggregate_samples(samples, method)` incl. empty handling
  - [X] `interpolate_empty_bins(power, empty_idx)` (incl. leading/trailing cases)
- [X] Update `tests/synthetic/generate_simple_hg_session.m`
  - [X] Add `stim.t`
  - [X] Case `fs > 1/dt` (downsample)
  - [X] Case `fs < 1/dt` (upsample warning)

**Tests to write**
- [X] `tests/unit/test_bin_high_gamma.m`
  - [X] Exact bin match (`fs = 1/dt`)
  - [X] Downsampling for each method (`mean/median/rms/max`)
  - [X] Upsampling triggers warning
  - [X] Empty bins interpolated correctly (incl. edges)
  - [X] Single bin & trailing samples edge cases
  - [X] Output length == `length(stim.t)`
  - [X] Non-negativity preserved
  - [X] Methods produce distinct expected statistics
- [X] `scripts/dev/demo_step2_binning.m`
  - [X] Build multiple `dt`
  - [X] 2×2 subplot: Original, Mean, Median, RMS
  - [X] Print stats & warnings
  - [X] Print "✓ Step 2 Complete: Binning working"

**Success criteria**
- [X] Unit tests pass
- [X] Demo illustrates downsampling & method differences
- [X] Empty bins handled with warnings
- [X] Output aligns exactly to `stim.t`

---

## Phase 1 — Core Model & First Integration (Steps 3–6)

### Step 3 — Gaussian Likelihood (Identity Link)
**What to implement**
- [X] `src/model/neglogli_gaussian.m`
  - [X] Signature: `[nll, grad] = neglogli_gaussian(w, X, y, lambda, D)` (identity)
  - [X] Compute residuals `r = y - X*w`
  - [X] `nll = (1/(2n))*sum(r.^2) + (λ/2)*sum((D*w).^2)`
  - [X] `grad = -(1/n)*X'*r + λ*D'*D*w` (use sparse ops; don't form `D'*D` explicitly)
  - [X] Normalize by `n`
  - [X] `validate_glm_inputs(...)` helper for dims/finite/λ≥0
- [X] No dependencies on Steps 1–2

**Tests to write**
- [X] `tests/unit/test_neglogli_gaussian.m`
  - [X] λ=0 matches manual SSE; gradient matches numerical
  - [X] λ>0 raises NLL; gradient includes penalty
  - [X] Finite-diff gradient checks at multiple `w`
  - [X] Edge cases: all zeros; perfect fit; dim mismatch errors
  - [X] Sparse vs dense X produce same results; large sparse memory sanity

**Demo**
- [X] `scripts/dev/demo_step3_likelihood.m`
  - [X] Small synthetic X,y, true `w`
  - [X] Compare NLL at true/zero/random `w`
  - [X] Plot NLL vs dims; overlay gradient arrows
  - [X] Finite-diff vs analytic gradient (max rel err < 1e-4)
  - [X] Print "✓ Step 3 Complete: Likelihood functions working"

**Success criteria**
- [X] Tests pass with high-precision gradient match
- [X] Penalty increases NLL
- [X] Sparse pathway correct
- [X] Validation catches bad inputs

---

### Step 4 — Basic Fitting Infrastructure
**What to implement**
- [ ] `src/model/fit_glm_map_hg.m`
  - [X] Signature: `fit = fit_glm_map_hg(X, y, cfg)`
  - [X] Init `w0=zeros(p,1)`
  - [X] fminunc (quasi-newton/trust-region) with gradient on
  - [X] Use `cfg.model.lambda`, `cfg.model.penalty_operator=D`
  - [X] Link: `'identity'` (this step)
  - [X] Options: `max_iter`, `tol_fun`, `tol_grad`, display
  - [X] Return: `w`, `mu = X*w`, `nll`, `exitflag`, `output`, `converged`, `grad_norm`
- [X] `src/model/predict_power.m`
  - [X] Identity: `mu = X*w`
  - [X] Log (stub for future): `mu = exp(X*w)`; clip negatives w/ warning for identity
- [X] Convergence checker
  - [X] `check_convergence(fit, cfg)` with warnings: not converged, large grad, NaN/Inf NLL

**Tests to write**
- [X] `tests/unit/test_fit_glm_map_hg.m`
  - [X] Simple 2D problem recovers `w_true`
  - [X] λ=0 matches least squares (`X\y`)
  - [X] λ>0 shrinks weights; greater λ → more shrinkage
  - [X] Perfect data converges quickly with low NLL
  - [X] Diagnostics captured: `exitflag`, `grad_norm`; warn on non-convergence
- [X] `tests/synthetic/generate_simple_fit_problem.m` (n=100, p=5, SNR≈3)

**Demo**
- [X] `scripts/dev/demo_step4_fitting.m`
  - [X] n=200, 3 features, known `w_true`, add noise
  - [X] Fit λ∈{0, 0.1, 1.0}; print `w_true` vs fitted
  - [X] Convergence status & metrics
  - [X] 2-panel plot: shrinkage curve; obs vs pred scatter (best λ)
  - [X] Print iterations, `grad_norm`, final NLL
  - [X] Print “✓ Step 4 Complete: Fitting infrastructure working”

**Success criteria**
- [X] Tests pass; accurate recovery (≤10% error at low noise)
- [X] Clear λ-shrinkage behavior
- [X] Convergence diagnostics populated
- [X] Pred vs obs corr > 0.9 (simple synthetic)

---

### Step 5 — Design Matrix Assembly (HG Wrapper)
**What to implement**
- [X] `src/features/assemble_design_matrix_hg.m`
  - [X] Set `cfg_hg = cfg`; set `include_spike_history=false`, clear history window
  - [X] Call existing `assemble_design_matrix(streams, states, sps, cfg_hg, stim)`
  - [X] `y = hg_binned.power`
  - [X] Validate `size(X,1)==numel(y)`; no NaN/Inf in X or y
  - [X] Return `X`, `y`, `colmap` (no spike history fields)
- [X] `validate_design_matrix_hg(X,y,colmap)`
  - [X] Confirm dims; verify no spike history cols
  - [X] Reasonable sparsity; expected regressor types present
- [X] `tests/synthetic/generate_single_event_hg.m`
  - [X] Single perceived call, known Gaussian bump kernel on flat baseline

**Tests to write**
- [X] `tests/unit/test_assemble_design_matrix_hg.m`
  - [X] Single perceived event → heard columns with causal window [0,2s]
  - [X] Single produced event → produced columns with symmetric window [−2s, +3s]
  - [X] No events → intercept-only, no error
  - [X] `colmap` fields & indices correct; dims consistent
- [X] `tests/e2e/test_design_matrix_integration.m`
  - [X] Assemble DM, fit with Step 4, fit completes; basic weight sanity

**Demo**
- [X] `scripts/dev/demo_step5_design_matrix.m`
  - [X] Build single perceived session (baseline=5, call at 2s, peak +2, σ=0.3s, delay 0.1s, dur 6s)
  - [X] `spy` plot of `X`; list regressors from `colmap`
  - [X] Show sample rows pre/during/post event
  - [X] Print "✓ Step 5 Complete: Design matrix assembly working"

**Success criteria**
- [X] Tests pass; correct dims and structure
- [X] No spike history columns
- [X] End-to-end fitting works on this DM

---

### Step 6 — First End-to-End Integration
**What to implement**
- [ ] Upgrade `tests/synthetic/generate_single_event_hg.m` (robust)
  - [ ] Duration 10s; dt=0.01; baseline=10
  - [ ] Single perceived at 5s
  - [ ] True kernel: Gaussian (peak=5, σ=0.2s, delay=0.15s)
  - [ ] Noise std=0.5
  - [ ] Save `hg_data`, `events`, `stim`, `cfg`, `ground_truth` to `tests/data/synthetic_single_event.mat`
- [ ] `src/model/extract_kernels.m`
  - [ ] Parse `w` using `colmap`; return struct per regressor with time vectors

**Tests to write**
- [ ] `tests/unit/test_extract_kernels.m`
  - [ ] Known `w`/`colmap` → correct parsing & time vectors
  - [ ] Missing regressors handled
- [ ] `tests/e2e/test_single_event_recovery.m`
  - [ ] Run full pipeline
  - [ ] Assert: converged; kernel vs truth corr > 0.8; peak bin ±1; peak mag within 30%; overall R² > 0.5

**Demo**
- [ ] `scripts/dev/demo_step6_first_integration.m`
  - [ ] Full pipeline with verbose banners
  - [ ] 2×2 figure: obs vs pred (with events & R²); scatter; fitted vs true kernel; residuals over time
  - [ ] Print convergence & kernel recovery metrics
  - [ ] Print “✓ CHECKPOINT 1: First end-to-end integration successful!”

**Success criteria**
- [ ] E2E tests pass with metrics thresholds
- [ ] Demo visuals match expectations
- [ ] No convergence warnings

---

## Phase 2 — Multiple Events (Step 7)

### Step 7 — Multi-Event Handling
**What to implement**
- [ ] `tests/synthetic/generate_multi_event_hg.m`
  - [ ] Duration 20s; dt=0.01; baseline=10; noise std=0.5
  - [ ] Events (non-overlapping): heard_addressed @ [3,9,15]; heard_overheard @ [5,12]; produced @ [7,17]
  - [ ] True kernels:
    - [ ] addressed: fast positive (peak=4 @0.1s, width=0.15s)
    - [ ] overheard: slow positive (peak=2 @0.3s, width=0.4s)
    - [ ] produced: biphasic (dip=−1 @−0.1s; peak=3 @0.2s)
  - [ ] Add ground truth struct; save to `tests/data/synthetic_multi_event.mat`
- [ ] Ensure `extract_kernels` returns separate fields for all three types
- [ ] `src/eval/compare_kernels.m` (metrics table per kernel)
  - [ ] Corr, peak error, magnitude ratio; pretty print table

**Tests to write**
- [ ] `tests/e2e/test_multi_event_recovery.m`
  - [ ] Fit with λ≈0.05
  - [ ] For each kernel: corr > 0.7; peak timing within 2 bins; peak mag within 40%
  - [ ] Overall R² > 0.6

**Demo**
- [ ] `scripts/dev/demo_step7_multi_event.m`
  - [ ] Event count summary
  - [ ] 2×2 figure: full trace w/ markers & R²; overlaid kernels (fitted vs true); bar plot of correlations; residual analysis
  - [ ] Kernel metrics table printed
  - [ ] Print “✓ CHECKPOINT 2: Multi-event kernel recovery successful!”

**Success criteria**
- [ ] All kernel thresholds met
- [ ] R² > 0.6
- [ ] Distinct, non-interfering kernels

---

## Phase 3 — State Effects (Step 8)

### Step 8 — Conversational State Effects
**What to implement**
- [ ] `tests/synthetic/generate_state_effect_hg.m`
  - [ ] 30s with epochs: [0–10]=spon (baseline=8), [10–20]=convo (baseline+4), [20–30]=spon
  - [ ] 2 perceived calls per epoch; kernel: Gaussian (peak=3, delay=0.15s, width=0.2s)
  - [ ] `state_convo=1` during [10–20], `state_spon=1` otherwise
  - [ ] True coefficients: β_convo=+4, β_spon=0; add ground truth
- [ ] Verify DM includes state columns (scalar indicators)
- [ ] `extract_kernels` returns `state_coeffs` map

**Tests to write**
- [ ] `tests/e2e/test_state_recovery.m`
  - [ ] λ≈0.05; kernel corr > 0.7
  - [ ] β_convo within 20% of +4; β_spon ~0 (±0.5)
  - [ ] R² > 0.7; mean predicted power convo > spon
- [ ] `tests/unit/test_state_in_design_matrix.m`
  - [ ] State columns are 0/1, mutually exclusive
  - [ ] `colmap` correctly identifies state columns

**Demo**
- [ ] `scripts/dev/demo_step8_states.m`
  - [ ] State stats (%time each state); DM check
  - [ ] 2×3 figure: traces by epoch; kernel fit; state coefficients; baseline by state
  - [ ] State estimation table
  - [ ] Print “✓ CHECKPOINT 3: State effects working!”

**Success criteria**
- [ ] State coefficients accurate within tolerance
- [ ] R² improvement vs model without states
- [ ] Kernel consistent across states

---

## Phase 4 — Basis & Smoothness (Step 9)

### Step 9 — Basis Projection and Smoothness Regularization
**What to implement**
- [ ] `tests/synthetic/generate_noisy_kernel_hg.m`
  - [ ] 15s; 5 perceived; baseline=10; noise std=1.5 (SNR≈1)
  - [ ] True smooth Gaussian kernel (peak=4, delay 0.15s, width 0.25s)
  - [ ] Version A: fine sampling (raw lags)
  - [ ] Version B: 10 raised-cosine basis funcs spanning window
  - [ ] Save ground truth
- [ ] Verify `build_basis_block.m` integration with HG wrapper
- [ ] Verify `smoothness_penalty.m` D matrix with λ sweep
- [ ] `src/eval/evaluate_kernel_smoothness.m`
  - [ ] Roughness = sum(diff(kernel,2).^2); MSE to truth; return trade-offs

**Tests to write**
- [ ] `tests/unit/test_basis_smoothness.m`
  - [ ] Basis projection reconstructs smooth kernels with small error
  - [ ] λ∈[0,0.01,0.1,1.0]: roughness decreases with λ; MSE shows U-shape
- [ ] `tests/e2e/test_basis_recovery.m`
  - [ ] Compare configs: no basis/no pen; basis/no pen; basis+λ=0.1
  - [ ] Assert basis+pen has best MSE to truth; report corr/MSE/roughness

**Demo**
- [ ] `scripts/dev/demo_step9_basis_smoothness.m`
  - [ ] Fit 4 configs: no basis/λ=0; basis/λ=0; basis/λ=0.01; basis/λ=0.1
  - [ ] 2×2 figure: kernels; MSE vs log λ (U-shape); roughness; traces & R²
  - [ ] Comparison table printed; flag best
  - [ ] Print “✓ CHECKPOINT 4: Regularization working!”

**Success criteria**
- [ ] Param count reduced via basis
- [ ] Roughness decreases with λ
- [ ] Optimal λ interior; regularized kernel best MSE & interpretable

---

## Phase 5 — Log Link (Step 10)

### Step 10 — Log Link Function
**What to implement**
- [ ] Update `src/model/neglogli_gaussian.m`
  - [ ] Add `'link'` param (default `'identity'`)
  - [ ] If `'log'`: `η=X*w`, `μ=exp(η)` (clip `η<=20`), same NLL form, gradient `-(1/n)X'((y-μ).*μ) + λ D' D w`
  - [ ] Guard NaN/Inf in `μ`
- [ ] Update `src/model/predict_power.m`
  - [ ] Implement `μ=exp(X*w)` for `'log'` with safety clip
- [ ] Update `src/model/fit_glm_map_hg.m`
  - [ ] Pass link to likelihood/predict; validate link ∈ {`identity`,`log`}
- [ ] `tests/synthetic/generate_multiplicative_hg.m`
  - [ ] Multiplicative truth: `log μ = log(baseline) + β_h * k_h + β_p * k_p`
  - [ ] Duration 20s; 4 heard, 3 produced; some co-occurrence; save truth & link

**Tests to write**
- [ ] `tests/unit/test_log_link.m`
  - [ ] Identity vs log NLL on same data (log uses exp correctly)
  - [ ] Finite-diff gradient matches
  - [ ] Stability for large positive `w` (no overflow; clipping works)
  - [ ] Log predictions strictly positive
- [ ] `tests/e2e/test_log_link.m`
  - [ ] Fit both links on multiplicative session
  - [ ] Log link higher R²; more homoscedastic residuals; better co-occurrence modeling

**Demo**
- [ ] `scripts/dev/demo_step10_log_link.m`
  - [ ] 2×3 figure: fits & R²; residual diagnostics; co-occurrence window (obs vs identity vs log)
  - [ ] Comparison metrics table; usage guidance
  - [ ] Print “✓ CHECKPOINT 5: Log link working!”

**Success criteria**
- [ ] Likelihood/gradient correct & stable
- [ ] Log link superior on multiplicative data; predictions positive

---

## Phase 6 — Cross-Validation (Step 11)

### Step 11 — Blocked K-Fold Cross-Validation
**What to implement**
- [ ] `src/model/crossval_blocked_hg.m`
  - [ ] Block K contiguous folds; λ grid (log-spaced) from cfg
  - [ ] For each λ, each fold: train on K−1, test on held-out
  - [ ] Metric: MSE or NLL; average & SE across folds
  - [ ] Return: `lambda_grid`, `test_error`, `test_error_se`, `best_lambda`, `best_idx`, `fold_errors`, `convergence_flags`
- [ ] `src/utils/generate_blocked_folds.m`
  - [ ] Split 1…n into K contiguous blocks (handle n%K≠0)
- [ ] `src/plot/plot_cv_curve.m`
  - [ ] Plot test error vs log λ, error bars, mark optimum

**Tests to write**
- [ ] `tests/unit/test_crossval_hg.m`
  - [ ] Fold coverage & contiguity (n=100,K=5)
  - [ ] λ grid count & log spacing correct
  - [ ] U-shaped error curve; best λ interior
  - [ ] Edge cases: small n; large K; single λ
- [ ] `tests/e2e/test_cv_integration.m`
  - [ ] Session with known optimal λ; select near truth; kernels recovered; U-shape curve

**Demo**
- [ ] `scripts/dev/demo_step11_crossval.m`
  - [ ] K=5; λ grid 20 values [0.001,10]; metric MSE
  - [ ] Progress prints per λ and fold
  - [ ] 2×2 figure: CV curve; zoom near optimum; per-fold curves; fit at optimal λ
  - [ ] Printed summary & model comparison
  - [ ] Print “✓ CHECKPOINT 6: Cross-validation working!”

**Success criteria**
- [ ] CV completes; best λ interior
- [ ] Reasonable fold variability
- [ ] Fit at best λ yields good recovery

---

## Phase 7 — Robustness (Step 12)

### Step 12 — Robustness & Edge Cases
**What to implement**
- [ ] `src/utils/validate_config_hg.m`
  - [ ] Check: dt>0; kernel windows valid; optim tolerances/iters; CV settings; link ∈{identity,log}; binning method valid
  - [ ] Apply defaults; throw informative errors
- [ ] Enhance validations/warnings across code
  - [ ] `load_high_gamma.m`: missing fields; malformed `channel_ids`; `fs>0`
  - [ ] `bin_high_gamma.m`: empty windows; warn if >10% bins empty; error if `hg.t` doesn’t overlap `stim.t`
  - [ ] `fit_glm_map_hg.m`: NaN/Inf in X or y; high condition number warning; actionable convergence suggestions
  - [ ] `crossval_blocked_hg.m`: error if K>n; warn if fold size<50; handle failed folds
- [ ] Comprehensive edge tests in one suite

**Tests to write**
- [ ] `tests/unit/test_edge_cases_hg.m`
  - [ ] Empty events → intercept-only, succeeds
  - [ ] Single event minimal case → warning but succeeds
  - [ ] Events at boundaries handled
  - [ ] Very short session; very long session (memory)
  - [ ] High noise (SNR < 0.5) converges with poor R² (expected)
  - [ ] Perfectly correlated regressors → collinearity warning
  - [ ] All-zero power → error
  - [ ] NaNs in data → error
  - [ ] Negative power → clipped with warning
- [ ] `tests/integration/test_error_messages.m`
  - [ ] Each error type prints what/where/how to fix (e.g., λ boundary → suggest expanded grid)

**Demo**
- [ ] `scripts/dev/demo_step12_robustness.m`
  - [ ] Run scenario suite: Empty, Minimal, High noise, Collinearity, Boundary λ, Invalid inputs
  - [ ] Diagnostic summary table printed
  - [ ] Print “✓ CHECKPOINT 7: Robustness validated! Production-ready.”

**Success criteria**
- [ ] All edge tests pass; informative errors/warnings
- [ ] No silent failures; early validation

---

## Phase 8 — Full Pipeline (Step 13)

### Step 13 — Full Pipeline Integration Script
**What to implement**
- [ ] `scripts/run_fit_high_gamma.m` (main entry)
  - [ ] Stage 1: Setup & Validation (config, output dir, logging)
  - [ ] Stage 2: Data Loading (HG, events; compatibility checks; summary)
  - [ ] Stage 3: Preprocessing (time base, bin HG, streams, states; summary)
  - [ ] Stage 4: Design Matrix (assemble; basis optional; penalties; visualize; optional save)
  - [ ] Stage 5: Model Fitting (CV optional; select λ; fit final; extract kernels)
  - [ ] Stage 6: Evaluation (predictions; metrics; diagnostics; optional permutation)
  - [ ] Stage 7: Output Generation (`fit_results.mat`, plots, `summary_report.txt`)
  - [ ] Progress messages, timing, error handling, intermediate saves
- [ ] `tests/synthetic/generate_realistic_hg_session.m`
  - [ ] 120s; realistic conversational bouts; events: 15 addressed, 8 overheard, 12 produced (some clustering)
  - [ ] Baseline=10; state +3 convo; noise std=1.0
  - [ ] Plausible kernels; save as `tests/data/synthetic_realistic.mat`
- [ ] Output organization
  - [ ] Results dir `results/high_gamma_<session>_<timestamp>/`
  - [ ] Save config, `fit_results.mat`, figures (DM, CV, kernels, fit quality, diagnostics), summary report

**Tests to write**
- [ ] `tests/e2e/test_full_pipeline.m`
  - [ ] Run via `run_fit_high_gamma` on realistic synthetic
  - [ ] Files exist; fields present; plots generated
  - [ ] Metrics reasonable: R² > 0.4; kernel corr > 0.6; converged; CV interior λ

**Demo**
- [ ] `scripts/dev/demo_step13_full_pipeline.m`
  - [ ] Full verbose run; printed stage checklist & metrics
  - [ ] Open figures; show output directory tree
  - [ ] Print “✓ CHECKPOINT 8: Full pipeline working! Ready for real data.”

**Success criteria**
- [ ] End-to-end run stable; outputs complete & interpretable
- [ ] Performance targets met; processing time reasonable

---

## Phase 9 — Statistical Inference (Step 14)

### Step 14 — Permutation Testing
**What to implement**
- [ ] `src/model/perform_permutation_test_hg.m`
  - [ ] For `n_perms`: shuffle events, rebuild DM, refit with same λ, store kernels
  - [ ] Build null distributions per kernel timepoint
  - [ ] Compute 95% CIs; p-values; significant bins
  - [ ] Multiple-comparison correction (Bonferroni or FDR)
  - [ ] Return: `null_kernels`, `confidence_intervals`, `p_values`, `significant_bins`, `correction_method`
- [ ] Utilities
  - [ ] `src/utils/shuffle_event_labels.m` (preserve counts; break temporal structure)
  - [ ] `src/utils/compute_fdr_threshold.m` (Benjamini–Hochberg)
- [ ] Visualization
  - [ ] `src/plot/plot_kernels_with_confidence.m` (shaded CIs; significance marks; null overlays)

**Tests to write**
- [ ] `tests/unit/test_permutation_hg.m`
  - [ ] Shuffling preserves counts; breaks structure; randomness verified
  - [ ] Null distribution construction (mean≈0 for centered; std>0)
  - [ ] Strong effect → significant; pure noise → not significant
- [ ] `tests/e2e/test_permutation_hg.m`
  - [ ] Realistic session; `n_perms≈200`
  - [ ] Permutation completes; CIs narrower than null range
  - [ ] Some bins significant; p in [0,1]; correction applied

**Demo**
- [ ] `scripts/dev/demo_step14_permutation.m`
  - [ ] Progress bar for `n_perms`
  - [ ] 2×3 figure: kernels with CIs & significance; null peak hist; p-value heatmap; % significant bins
  - [ ] Statistical summary table printed
  - [ ] Methodological note printed
  - [ ] Print “✓ CHECKPOINT 9: Permutation testing complete! ✓ PROJECT COMPLETE”

**Success criteria**
- [ ] Runtime reasonable (~200 perms)
- [ ] Null centered near zero; observed outside null
- [ ] Correct p-values & multiple-comparison control
- [ ] Clear, publication-quality plots

---

## Global (from Spec) — Deliverables & Readiness

**Code Modules**
- [ ] All new HG functions in `src/` implemented
- [ ] Main pipeline script (`scripts/run_fit_high_gamma.m`)
- [ ] Demo scripts for all phases under `scripts/dev/`

**Testing**
- [ ] Unit tests for all core functions
- [ ] E2E tests for pipeline validation
- [ ] Synthetic data generators
- [ ] `./scripts/run_matlab_tests.sh` passes all tests

**Documentation**
- [ ] Spec checked in
- [ ] Implementation prompts (14 steps) checked in
- [ ] `config/defaults_hg.json` with comments
- [ ] README usage examples

**Validation Checkpoints**
- [ ] CP1: I/O & preprocessing working (after Step 2)
- [ ] CP2: Core model working; first kernel recovered (after Step 6)
- [ ] CP3: Multi-event recovery (after Step 7)
- [ ] CP4: Regularization validated (after Step 9)
- [ ] CP5: Log link working (after Step 10)
- [ ] CP6: CV working (after Step 11)
- [ ] CP7: Robustness validated (after Step 12)
- [ ] CP8: Full pipeline ready for real data (after Step 13)
- [ ] CP9: Permutation inference complete (after Step 14)

**Real Data Acceptance (MVP v0.1)**
- [ ] Pipeline runs on real session without errors
- [ ] R² > 0.05 on real data
- [ ] Kernels physiologically plausible
- [ ] Diagnostic plots reasonable
- [ ] CV selects interior λ
- [ ] All plots auto-generated

**Enhanced v0.2**
- [ ] Permutation test CIs included
- [ ] Gamma GLM alternative available
- [ ] Multiple sessions tested
- [ ] Documentation complete
