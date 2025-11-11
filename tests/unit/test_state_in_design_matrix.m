function tests = test_state_in_design_matrix()
%TEST_STATE_IN_DESIGN_MATRIX Unit tests for state columns in design matrix.
%   Verifies that state regressors are correctly represented as scalar
%   indicators and properly mapped in colmap.

tests = functiontests(localfunctions);
end


function test_state_columns_are_binary(testCase)
%TEST_STATE_COLUMNS_ARE_BINARY Verify state columns contain only 0s and 1s.

% generate session with states
data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% extract state columns
testCase.verifyTrue(isfield(colmap, 'states'), 'Should have states in colmap');

state_cols = colmap.states.cols;
testCase.verifyEqual(length(state_cols), 2, 'Should have 2 state columns');

% check each state column
for i = 1:length(state_cols)
    col_idx = state_cols(i);
    col_values = full(X(:, col_idx));

    % values should be 0 or 1
    unique_vals = unique(col_values);
    testCase.verifyTrue(all(ismember(unique_vals, [0, 1])), ...
        sprintf('State column %d should contain only 0 and 1', i));
end

end


function test_states_are_mutually_exclusive(testCase)
%TEST_STATES_ARE_MUTUALLY_EXCLUSIVE Verify states don't overlap.
%   At each time point, exactly one state should be active.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

state_cols = colmap.states.cols;

% extract state columns
state_matrix = full(X(:, state_cols));

% sum across states at each time point should be 1
state_sum = sum(state_matrix, 2);
testCase.verifyTrue(all(state_sum == 1), ...
    'Exactly one state should be active at each time point');

end


function test_colmap_identifies_states(testCase)
%TEST_COLMAP_IDENTIFIES_STATES Verify colmap correctly identifies state columns.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% verify colmap structure
testCase.verifyTrue(isfield(colmap, 'states'), 'colmap should have states field');
testCase.verifyTrue(isfield(colmap.states, 'cols'), 'states should have cols field');
testCase.verifyTrue(isfield(colmap.states, 'names'), 'states should have names field');

% verify names
testCase.verifyEqual(length(colmap.states.names), 2, 'Should have 2 state names');
testCase.verifyTrue(ismember('convo', colmap.states.names), ...
    'Should have convo state');
testCase.verifyTrue(ismember('spon', colmap.states.names), ...
    'Should have spon state');

end


function test_state_columns_coverage(testCase)
%TEST_STATE_COLUMNS_COVERAGE Verify states cover expected time periods.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% extract state columns
state_names = colmap.states.names;
convo_idx = find(strcmp(state_names, 'convo'));
spon_idx = find(strcmp(state_names, 'spon'));

convo_col = colmap.states.cols(convo_idx);
spon_col = colmap.states.cols(spon_idx);

convo_active = full(X(:, convo_col));
spon_active = full(X(:, spon_col));

% conversational state should be active for ~1/3 of time [10-20s out of 30s]
pct_convo = 100 * sum(convo_active) / length(convo_active);
testCase.verifyGreaterThan(pct_convo, 30, 'Convo should be active >30% of time');
testCase.verifyLessThan(pct_convo, 40, 'Convo should be active <40% of time');

% spontaneous state should be active for ~2/3 of time
pct_spon = 100 * sum(spon_active) / length(spon_active);
testCase.verifyGreaterThan(pct_spon, 60, 'Spon should be active >60% of time');
testCase.verifyLessThan(pct_spon, 70, 'Spon should be active <70% of time');

end


function test_state_alignment_with_ground_truth(testCase)
%TEST_STATE_ALIGNMENT_WITH_GROUND_TRUTH Verify state columns match ground truth.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% extract state columns
state_names = colmap.states.names;
convo_idx = find(strcmp(state_names, 'convo'));
spon_idx = find(strcmp(state_names, 'spon'));

convo_col = colmap.states.cols(convo_idx);
spon_col = colmap.states.cols(spon_idx);

convo_active = logical(full(X(:, convo_col)));
spon_active = logical(full(X(:, spon_col)));

% compare to ground truth
gt_convo = data.ground_truth.convo_mask;
gt_spon = data.ground_truth.spon_mask;

testCase.verifyEqual(convo_active, gt_convo, ...
    'Convo state should match ground truth');
testCase.verifyEqual(spon_active, gt_spon, ...
    'Spon state should match ground truth');

end


function test_states_are_scalar_not_kernels(testCase)
%TEST_STATES_ARE_SCALAR_NOT_KERNELS Verify states are not convolved kernels.
%   State columns should be simple indicators, not lag structures.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% states should have exactly 2 columns (convo, spon)
state_cols = colmap.states.cols;
testCase.verifyEqual(length(state_cols), 2, ...
    'States should be 2 columns, not lag structures');

% verify states don't have lag structure like kernels
testCase.verifyFalse(isfield(colmap.states, 'lag_bins'), ...
    'States should not have lag_bins like kernels');
testCase.verifyFalse(isfield(colmap.states, 'window_s'), ...
    'States should not have window_s like kernels');

end


function test_state_extraction(testCase)
%TEST_STATE_EXTRACTION Verify extract_kernels properly extracts state coeffs.

data = generate_state_effect_hg();
hg_binned = bin_high_gamma(data.hg_data, data.stim, data.cfg);
[X, y, colmap] = assemble_design_matrix_hg(hg_binned, data.events, data.stim, data.cfg);

% fit model
fit = fit_glm_map_hg(X, y, data.cfg);

% extract kernels and states
kernels = extract_kernels(fit.w, colmap, data.stim);

% verify states extracted
testCase.verifyTrue(isfield(kernels, 'states'), ...
    'extract_kernels should extract states');
testCase.verifyTrue(isfield(kernels.states, 'convo'), ...
    'Should have convo state coefficient');
testCase.verifyTrue(isfield(kernels.states, 'spon'), ...
    'Should have spon state coefficient');

% verify state coefficients are scalars
testCase.verifyTrue(isscalar(kernels.states.convo), ...
    'Convo coefficient should be scalar');
testCase.verifyTrue(isscalar(kernels.states.spon), ...
    'Spon coefficient should be scalar');

end
