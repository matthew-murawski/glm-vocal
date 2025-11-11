GLM for High Gamma Power in Natural Marmoset Conversation — Developer Specification
Status: v0.2 (Revised with dual verification approach)Owner: Zhao LabLanguage: MATLABGoal: Adapt the existing spike-based GLM pipeline to predict continuous high gamma power traces during naturalistic marmoset conversations, reusing preprocessing and feature engineering infrastructure while replacing the Poisson model with a Gaussian (or Gamma) GLM.

0) Executive Summary (TL;DR)
What's the same (reuse ~85% of existing pipeline):
* All preprocessing: time base construction, event stream generation, conversational state detection
* All feature engineering: lagged kernels, raised-cosine basis projection, design matrix assembly, smoothness penalties
* All cross-validation logic: blocked CV, fold generation, λ selection
* All visualization and I/O: kernel plots, design matrix previews, config handling
What's different (core statistical model):
* Input data: Continuous high gamma power trace [nT x 1] instead of spike counts
* Distribution family: Gaussian (or Gamma) instead of Poisson
* Link function: Identity (linear) or log, instead of log-only
* Likelihood: Sum of squared errors (Gaussian-identity) or Gaussian/Gamma log-likelihood
* Outputs: Predicted power trace with R², MSE, correlation metrics instead of Poisson deviance
Implementation strategy:
Test-Driven Development (TDD) with synthetic data first, building in modular chunks so you see each component working before integration. Start with trivial synthetic sessions, gradually add complexity, then move to real data.
Dual verification approach:
* Machine validation: Unit and e2e tests using MATLAB test suite, runnable via ./scripts/run_matlab_tests.sh
* Human validation: Demo scripts in scripts/dev/ with verbose output and visualizations for manual inspection at each development stage
See accompanying HIGH_GAMMA_GLM_IMPLEMENTATION_PROMPTS.md for detailed 14-step implementation plan with explicit checkpoints.

1) Scope & Requirements
1.1 Current Scope
Single-session, single-channel high gamma power trace prediction using the same conversational event regressors as the spike-based GLM.
Inputs:
* High gamma power trace (continuous, positive values, one sample per time bin)
* Labeled produced calls (onset, offset, optional type)
* Labeled perceived calls (onset, offset, optional type)
* Session metadata (subject, date, recording parameters)
Regressors (identical to spike GLM):
* heard_addressed, heard_overheard, heard_any: binary streams with causal kernels [0, 2s]
* produced_*: split by context (spontaneous/after_heard/after_produced) or call_type, symmetric kernels [-2s, +3s]
* state_convo, state_spon: scalar context flags
* No spike history (not applicable for continuous power)
Model:
* Gaussian GLM with identity or log link (user-selectable)
* Alternative: Gamma GLM with log link for strictly positive skewed data
* L2 smoothness penalties on kernel blocks (identical to spike GLM)
* Blocked K-fold cross-validation to select λ
Outputs:
* Fitted kernels showing temporal dynamics of high gamma response to calls
* Predicted high gamma power trace
* R², MSE, correlation, explained variance
* Permutation test confidence intervals for kernels
* Standard plots (kernels, predicted vs observed, CV curve)
1.2 Out of Scope (for now)
* Multi-channel analysis (though architecture should allow future extension)
* Spectral decomposition (assume high gamma power is pre-computed from LFP)
* Non-linear kernels or multiplicative interactions between regressors
* Hierarchical/mixed-effects across sessions
1.3 Key Differences from Spike GLM
Aspect	Spike GLM	High Gamma GLM
Response variable	Spike counts (non-negative integers)	High gamma power (continuous, positive)
Distribution	Poisson	Gaussian or Gamma
Link function	log (only)	identity or log (user choice)
Spike history	Yes (causal kernel)	No (not applicable)
Typical R²	0.05–0.30	Potentially higher (0.10–0.50+)
Preprocessing	Bin spike times → counts	Downsample/bin power trace
2) Repository Structure (extends existing glm-vocal/)
New Files for High Gamma
glm-vocal/
├─ src/
│  ├─ io/
│  │  └─ load_high_gamma.m              # NEW: load .h5/.mat high gamma trace
│  │  └─ validate_inputs_hg.m           # NEW: validation for high gamma data
│  ├─ preprocess/
│  │  └─ bin_high_gamma.m               # NEW: downsample/bin continuous trace to stim.t
│  ├─ features/
│  │  └─ assemble_design_matrix_hg.m    # NEW: wrapper that excludes spike history
│  ├─ model/
│  │  ├─ neglogli_gaussian.m            # NEW: Gaussian NLL + gradient
│  │  ├─ neglogli_gamma.m               # NEW: Gamma NLL + gradient (optional)
│  │  ├─ fit_glm_map_hg.m               # NEW: wrapper using Gaussian likelihood
│  │  ├─ predict_power.m                # NEW: prediction function, returns μ
│  │  ├─ crossval_blocked_hg.m          # NEW: CV with MSE metric
│  │  └─ perform_permutation_test_hg.m  # NEW: statistical inference on kernels
│  ├─ eval/
│  │  └─ metrics_hg.m                   # NEW: R², MSE, correlation metrics
│  └─ plot/
│     └─ plot_power_diagnostics.m       # NEW: residual diagnostics
├─ scripts/
│  ├─ run_fit_high_gamma.m              # NEW: main entry point for HG analysis
│  └─ dev/                               # NEW: human-checkable demo scripts
│     ├─ demo_step1_io.m
│     ├─ demo_step2_binning.m
│     ├─ demo_step3_likelihood.m
│     ├─ demo_step4_fitting.m
│     ├─ demo_step5_design_matrix.m
│     ├─ demo_phase1_single_event.m      # CHECKPOINT: Phase 1 complete
│     ├─ demo_phase2_multi_event.m       # CHECKPOINT: Phase 2 complete
│     ├─ demo_phase3_states.m            # CHECKPOINT: Phase 3 complete
│     ├─ demo_phase4_basis_smoothness.m  # CHECKPOINT: Phase 4 complete
│     ├─ demo_phase5_log_link.m          # CHECKPOINT: Phase 5 complete
│     ├─ demo_phase6_crossval.m          # CHECKPOINT: Phase 6 complete
│     ├─ demo_phase7_robustness.m        # CHECKPOINT: Phase 7 complete
│     ├─ demo_phase8_full_pipeline.m     # CHECKPOINT: Phase 8 ready
│     └─ demo_phase9_permutation.m       # CHECKPOINT: Phase 9 complete, project done
└─ tests/
   ├─ synthetic/                         # NEW: synthetic data generators
   │  ├─ generate_simple_hg_session.m
   │  ├─ generate_single_event_hg.m
   │  ├─ generate_multi_event_hg.m
   │  ├─ generate_state_effect_hg.m
   │  ├─ generate_noisy_kernel_hg.m
   │  └─ generate_multiplicative_hg.m
   ├─ unit/                              # NEW: unit tests for HG components
   │  ├─ test_load_high_gamma.m
   │  ├─ test_bin_high_gamma.m
   │  ├─ test_neglogli_gaussian.m
   │  ├─ test_fit_glm_map_hg.m
   │  ├─ test_assemble_design_matrix_hg.m
   │  ├─ test_log_link.m
   │  ├─ test_crossval_hg.m
   │  ├─ test_basis_smoothness.m
   │  └─ test_edge_cases_hg.m
   └─ e2e/                               # NEW: end-to-end tests
      ├─ test_single_event_recovery.m
      ├─ test_multi_event_recovery.m
      ├─ test_state_recovery.m
      └─ test_permutation_hg.m
Strategy: Keep all existing spike GLM code intact. Create parallel _hg versions only for the statistical model layer. Reuse all preprocessing and feature code by calling existing functions.

3) Data Handling & Formats
3.1 High Gamma Power Trace
MATLAB .mat file produced by export_continuous_hg_for_glm.
* Expected variables and fields:
    * lfp_data.lfp_data — double matrix [nSamples × nChannels] of high gamma envelope (downsampled; non-negative).
    * lfp_data.lfp_fs — scalar sampling rate (Hz) of lfp_data.lfp_data.
    * lfp_data.lfp_t0 — scalar start time in audio coordinates (seconds) and/or
    * metadata.timestamps — double vector [nSamples × 1] of audio-time stamps for each sample (preferred if present).
    * lfp_data.channel_ids — vector [nChannels × 1] mapping columns to channel IDs.
    * (Optional) metadata struct — includes processing options, session info; used for provenance.
* load_high_gamma.m must:
    1. Load the MAT file;
    2. Select a single channel according to config (e.g., channel_id, or reducer strategy like median across channels);
    3. Return a struct:
        * power — [nSamples × 1] double (selected or reduced channel);
        * fs — scalar (use lfp_data.lfp_fs);
        * t — [nSamples × 1] double time vector. Prefer metadata.timestamps; otherwise construct t = lfp_data.lfp_t0 + (0:nSamples-1)'/lfp_data.lfp_fs.
        * session_id, channel_id (optional for logging).
3.2 Binning/Downsampling High Gamma to stim.t
High gamma traces may be sampled at TargetFs (e.g., 50–200 Hz). GLM bins are typically at dt (e.g., 0.01 s → 100 Hz).
* Aggregation methods (configurable): mean (recommended), median, RMS, max.
* Mapping:
    * If fs ≥ 1/dt, aggregate samples within each bin (no upsampling).
    * If fs < 1/dt, prefer keeping dt ≥ 1/fs in config; otherwise, allow linear interpolation to the bin centers (warn in logs).
* Handle edge cases: empty bins (interpolate), trailing samples, ensure output length matches stim.t.
* Output is non-negative.
3.3 Event Labels (unchanged from spike GLM)
Reuse existing load_labels.m and event processing. Format identical to spike GLM:
* Struct array with fields: kind ('produced'|'perceived'), t_on, t_off, label
* Support both MAT and Audacity TXT formats
* Timestamps in seconds relative to session start
3.4 Time Base (unchanged)
Reuse build_timebase.m. The stim.t grid is the same; only the response variable changes.
3.5 Regressor Streams (unchanged)
Reuse build_streams.m, compute_states.m, classify_produced_events.m. All feature engineering stays identical.

4) Feature Engineering → Design Matrix
4.1 Reuse Existing Infrastructure
No changes needed to core feature functions:
* build_kernel_block.m: creates lagged copies with causal/symmetric modes
* build_basis_block.m: raised-cosine basis projection for smooth kernels
* smoothness_penalty.m: second-difference penalty matrices for regularization
These functions are agnostic to whether predicting spikes or high gamma power.
4.2 Modified Assembly Function
Create thin wrapper assemble_design_matrix_hg.m that:
* Calls existing assemble_design_matrix.m
* Automatically excludes spike_history (not applicable for continuous power)
* Takes hg_binned instead of spike counts as response variable
* Returns same design matrix structure: X (sparse), y (hg_binned), colmap (column mapping)
Design matrix layout (column order):
[ intercept | heard_addressed | heard_overheard | heard_any | 
  produced_* (by context or type) | state_convo | state_spon ]

5) Statistical Models
5.1 Model Choice: Gaussian vs Gamma
Option 1: Gaussian GLM (recommended for MVP)
Conceptual model:
* Response follows normal distribution around predicted mean
* Mean can be linear (identity link) or exponential (log link) function of predictors
* Assumes additive Gaussian noise
When appropriate:
* High gamma power is approximately symmetric around mean
* No extreme outliers or heavy tails
* Variance roughly constant or stabilized by log link
Option 2: Gamma GLM (alternative for skewed data)
Conceptual model:
* Response follows Gamma distribution (positive, right-skewed)
* Variance scales with mean (realistic for power measurements)
* Typically uses log link
When appropriate:
* High gamma power is strictly positive and right-skewed
* Variance increases with mean power
* Occasional large values, long right tail
For MVP: Start with Gaussian. Add Gamma as enhancement if diagnostics show poor fit.
5.2 Gaussian GLM Mathematical Structure
Identity Link
Model: response = linear_combination_of_predictors + gaussian_noise
Objective function: sum_of_squared_residuals + penalty_on_kernel_smoothness
Gradient: linear in residuals, efficient to compute
Special case: Without penalty, this is ordinary least squares
Log Link
Model: response ~ Normal(exp(linear_combination_of_predictors), variance)
Purpose: Ensures predicted power is always positive
Gradient: Chain rule through exponential, involves both residuals and predictions
Numerical care: Clip linear predictor before exp to prevent overflow
5.3 Optimization Strategy
Approach:
* Use gradient-based optimizer (fminunc or L-BFGS)
* Initialize at zero (intercept-only model)
* Converge when objective and gradient change below tolerances
Numerical stability:
* Sparse matrix operations throughout
* Apply penalty as matrix-vector product, not explicit matrix
* For log link: clip values before exp

6) Cross-Validation & Model Selection
6.1 Blocked Cross-Validation
Purpose: Select optimal penalty strength λ
Procedure:
1. Split time into K contiguous folds (respects temporal structure)
2. For each candidate λ:
    * Train on K-1 folds, evaluate on held-out fold
    * Compute test metric (MSE or NLL)
3. Select λ minimizing average test metric
Key difference from spike GLM: Use MSE metric instead of Poisson deviance
6.2 Smoothness Penalty
Purpose: Regularize kernels to be temporally smooth
Implementation:
* Second-difference operator on kernel weights
* Penalizes curvature (discretized second derivative)
* Separate penalties for heard, produced kernels
* No penalty on intercept or state scalars

7) Outputs & Artifacts
7.1 Main Output File
Location: results/high_gamma_<session>_<timestamp>/fit_results.mat
Contents include:
* Configuration, input data, preprocessed streams
* Design matrix, penalty operators
* Cross-validation results, selected λ
* Fitted weights and unpacked kernels
* Predicted power trace
* Goodness-of-fit metrics
* Permutation test results (if performed)
7.2 Metrics
Computed for high gamma:
* R²: Fraction of variance explained
* MSE/RMSE: Mean squared/root mean squared error
* MAE: Mean absolute error
* Correlation: Between observed and predicted
* NLL per bin: Negative log-likelihood
7.3 Plots
Automatically generated:
1. Kernels: Temporal response profiles with confidence intervals
2. CV curve: Test metric vs penalty strength, shows selected λ
3. Fit quality: Observed vs predicted traces, scatter plot
4. Diagnostics: Residual histogram, Q-Q plot, heteroscedasticity check
5. Design matrix: Sparse structure visualization

8) Configuration
8.1 Configuration File
File: config/defaults_hg.json
Key sections:
* Time discretization: bin width (dt)
* Kernel windows: time ranges for each regressor type
* High gamma specific: binning method, link function, model type
* Cross-validation: number of folds, λ grid, evaluation metric
* Conversational state rules: response windows, gap thresholds
* Optimizer settings: tolerances, max iterations
* Basis projection: raised-cosine parameters

9) Test-Driven Development (TDD) Plan
9.1 Dual Verification Approach
Machine Validation:
* Unit tests for individual functions
* End-to-end tests for full pipeline
* Run via ./scripts/run_matlab_tests.sh
* Automated pass/fail criteria
* Fast execution for rapid iteration
Human Validation:
* Demo scripts in scripts/dev/ for each phase
* Verbose terminal output explaining steps
* Visualizations: plots of data, fits, kernels
* Clear success criteria printed with checkmarks
* Manual inspection builds intuition
Checkpoint Pattern:
1. Code implemented
2. Automated tests run and pass
3. Demo script created with verbose output
4. Developer runs demo, inspects plots
5. If good → proceed; if issues → debug
9.2 Synthetic Data Progression
Phase 0 (Steps 1-2): Foundation
* Flat high gamma trace (constant value)
* Load, bin, verify preservation
* Checkpoint: Can load and process basic data
Phase 1 (Steps 3-6): Single Event
* One perceived call, known kernel shape
* Fit GLM, recover kernel
* Checkpoint: First end-to-end integration working
Phase 2 (Step 7): Multi-Event
* Multiple events, different types
* Different kernels for addressed/overheard/produced
* Checkpoint: Handle complexity of real data structure
Phase 3 (Step 8): States
* Conversational vs spontaneous context effects
* Baseline differs by state
* Checkpoint: State coefficients work
Phase 4 (Step 9): Basis & Smoothness
* Noisy kernel, test regularization
* Checkpoint: Regularization improves recovery
Phase 5 (Step 10): Log Link
* Multiplicative effects
* Test link function choice
* Checkpoint: Both links work
Phase 6 (Step 11): Cross-Validation
* Automatic λ selection
* Checkpoint: Pipeline nearly complete
Phase 7 (Step 12): Robustness
* Edge cases, error handling
* Checkpoint: Production-ready validation
Phase 8 (Step 13): Full Pipeline
* Realistic synthetic session
* Complete end-to-end
* Checkpoint: Ready for real data
Phase 9 (Step 14): Permutation Testing
* Statistical inference
* Checkpoint: Project complete!
9.3 Major Checkpoints Summary
Checkpoint	After Step	Milestone	Success Indicator
1	Step 2	I/O working	Can load and bin high gamma
2	Step 4	Core model working	Simple fits converge
3	Step 6	Phase 1 complete	First kernel recovered
4	Step 9	Regularization validated	Basis+penalty improves fit
5	Step 11	Pipeline nearly done	CV working, all features present
6	Step 13	Ready for real data	Full pipeline on realistic synthetic
7	Step 14	Project complete	Permutation tests, production ready
10) Error Handling & Validation
10.1 Input Validation
Checks on high gamma data:
* Power is numeric vector, all finite
* Warn if negative values (set to zero)
* Error if all zeros or no variability
* Validate sampling rate
Error message format:
* Use MATLAB error identifiers
* Provide actionable guidance
10.2 Numerical Stability
In likelihood computation:
* Clip linear predictor before exp (log link)
* Check for overflow/underflow
In optimization:
* Monitor condition numbers
* Check gradient norm during iterations
* Use sparse matrix operations
10.3 Fit Diagnostics
Post-fit checks:
* Convergence flag
* Predictions are positive (log link)
* R² in reasonable range
* Residual diagnostics (histogram, Q-Q, autocorrelation)
10.4 Edge Case Handling
Common issues:
* Empty regressor streams → exclude from model
* Rare events → warn about low power
* Missing data → apply mask consistently
* Collinearity → warn, regularization helps
* Boundary λ selection → suggest wider grid

11) Deliverables Checklist
Code Modules
* [ ] All new high gamma functions in src/
* [ ] Main pipeline script
* [ ] Demo scripts for all phases in scripts/dev/
Testing
* [ ] Unit tests for all core functions
* [ ] E2E tests for pipeline validation
* [ ] Synthetic data generators
* [ ] All tests pass via ./scripts/run_matlab_tests.sh
Documentation
* [ ] This specification (conceptual framework)
* [ ] Implementation prompts (step-by-step instructions)
* [ ] Configuration file with comments
* [ ] README with usage examples
Validation Milestones
* [ ] Checkpoint 1: I/O and preprocessing working
* [ ] Checkpoint 2: Core model working
* [ ] Checkpoint 3: Phase 1 complete (first kernel recovery)
* [ ] Checkpoint 4: Regularization validated
* [ ] Checkpoint 5: CV working (pipeline nearly complete)
* [ ] Checkpoint 6: Ready for real data
* [ ] Checkpoint 7: Project complete (permutation testing)
Real Data Validation
* [ ] Pipeline runs on real session without errors
* [ ] R² > 0.05 (modest but non-zero)
* [ ] Kernels physiologically plausible
* [ ] Diagnostic plots show reasonable fit

12) Success Criteria
MVP Acceptance (v0.1):
* All synthetic tests pass (>90% success rate)
* Real data pipeline completes
* R² > 0.05 on real data
* Kernels smooth and interpretable
* Cross-validation selects interior λ
* All plots generate correctly
Enhanced Version (v0.2):
* Permutation test confidence intervals
* Gamma GLM alternative
* Tested on multiple sessions
* Complete documentation

13) Timeline Estimate
Development phases with estimated time:
Phase	Steps	Focus	Time
Phase 0	1-2	Foundation	3-4 hrs
Phase 1	3-6	Core model + integration	8-10 hrs
Phase 2-5	7-10	Add complexity	9-13 hrs
Phase 6-7	11-12	CV + robustness	6-8 hrs
Phase 8-9	13-14	Full pipeline + permutation	7-9 hrs
Total	14 steps	Complete	33-44 hrs
With debugging buffer: 40-53 hours total
Working schedules:
* Full-time: 1-1.5 weeks
* Part-time (20 hrs/week): 2-3 weeks
* Part-time (10 hrs/week): 4-5 weeks

14) References
Theoretical Background:
* McCullagh & Nelder (1989). Generalized Linear Models
* Friedman, Hastie, Tibshirani (2010). Regularization paths for GLMs
Implementation Guides:
* MATLAB GLM documentation
* Pillow et al. (2008). Raised-cosine basis reference
Existing Infrastructure:
* Your glm-vocal repository spike GLM implementation

Appendix A — Minimal Loader Contract for MAT Export (new)
src/io/load_high_gamma.m (MAT branch) — expected behavior:
Input: cfg.input.mat.path, cfg.input.mat.channel_strategy, cfg.input.mat.channel_id, cfg.input.mat.reduce.
Process:
1. S = load(path, 'lfp_data', 'metadata', 'lfp_fs', 'lfp_t0', 'channel_ids');
2. Determine channel vector C from lfp_data.lfp_data and lfp_data.channel_ids.
3. Select single channel (channel_id) or reduce across channels (median/mean).
4. Build time vector t from metadata.timestamps if present; otherwise t = lfp_data.lfp_t0 + (0:n-1)'/lfp_data.lfp_fs.
Output struct: hg.power, hg.fs, hg.t, hg.session_id (if available), hg.channel_id.
This completes the minimal, necessary adaptation to use export_continuous_hg_for_glm outputs without altering the rest of the GLM design.


END OF SPECIFICATION
This specification provides the complete conceptual framework for adapting your spike-based GLM to continuous high gamma power. All actual implementation details, step-by-step prompts, and code structures are in the accompanying HIGH_GAMMA_GLM_IMPLEMENTATION_PROMPTS.md document.