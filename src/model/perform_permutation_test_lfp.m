function ptest = perform_permutation_test_lfp(kernels, Xd, D, best_lambda, cfg, stim)
% section setup
% perform a permutation test to determine the statistical significance of lfp glm kernels.
% this function generates a null distribution of kernel peak amplitudes by repeatedly
% fitting the gaussian glm to data with shuffled lfp responses. it then compares the peak
% amplitude of the real kernels to this null distribution to calculate p-values.
% it also computes confidence intervals for the kernel weights.

% section initialize
% set the number of permutations and prepare storage for null distribution statistics.
if isfield(cfg, 'n_perms') && ~isempty(cfg.n_perms)
    n_perms = cfg.n_perms;
else
    n_perms = 250;
end

kernel_fields = fieldnames(kernels);
null_peaks = struct();
null_kernels = struct();

for i = 1:numel(kernel_fields)
    field = kernel_fields{i};
    if isstruct(kernels.(field)) && isfield(kernels.(field), 'weights')
        null_peaks.(field) = nan(n_perms, 1);
        null_kernels.(field) = nan(n_perms, numel(kernels.(field).weights));
    end
end

% section extract lfp parameters
% extract link function and optimizer config from cfg
if isfield(cfg, 'lfp') && isfield(cfg.lfp, 'link')
    link = cfg.lfp.link;
else
    link = 'identity';
end

if isfield(cfg, 'optimizer')
    optCfg = cfg.optimizer;
else
    optCfg = struct();
end

% section permutation loop
% iterate to build the null distribution.
fprintf('Running %d permutations for LFP null distribution...\n', n_perms);
for i = 1:n_perms
    if mod(i, 50) == 0
        fprintf('  Permutation %d/%d\n', i, n_perms);
    end

    % section shuffle data
    % circularly shift the lfp response to break temporal relationships.
    y_shuffled = circshift(Xd.y, randi(numel(Xd.y)));
    Xd_shuffled = Xd;
    Xd_shuffled.y = y_shuffled;

    % section fit model
    % fit the gaussian glm on the shuffled data to get a null kernel.
    wmap_null = fit_glm_gauss(Xd_shuffled, D, best_lambda, optCfg, link);
    kernels_null = unpack_params(wmap_null, Xd.colmap, cfg, stim);

    % section store results
    % extract and store the peak amplitude and weights of the null kernel.
    for j = 1:numel(kernel_fields)
        field = kernel_fields{j};
        if isstruct(kernels_null.(field)) && isfield(kernels_null.(field), 'weights')
            null_kernel_weights = kernels_null.(field).weights;
            null_peaks.(field)(i) = max(abs(null_kernel_weights));
            null_kernels.(field)(i, :) = null_kernel_weights;
        end
    end
end

% section calculate significance
% compute p-values and confidence intervals from the null distribution.
ptest = struct();
for i = 1:numel(kernel_fields)
    field = kernel_fields{i};
    if isstruct(kernels.(field)) && isfield(kernels.(field), 'weights')
        real_peak = max(abs(kernels.(field).weights));
        p_value = mean(null_peaks.(field) >= real_peak);

        ci_lower = prctile(null_kernels.(field), 2.5, 1);
        ci_upper = prctile(null_kernels.(field), 97.5, 1);

        ptest.(field) = struct('p_value', p_value, 'ci_lower', ci_lower, 'ci_upper', ci_upper);

        % report significance
        if p_value < 0.05
            fprintf('  %s: p=%.4f (significant)\n', field, p_value);
        else
            fprintf('  %s: p=%.4f (not significant)\n', field, p_value);
        end
    end
end

end
