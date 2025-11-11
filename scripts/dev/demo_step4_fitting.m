% demo_step4_fitting
% demo for step 4: basic fitting infrastructure for gaussian glm (identity)

% add paths if not already on path
if ~(exist('fit_glm_map_hg','file')==2)
    thisdir = fileparts(mfilename('fullpath'));
    reporoot = fullfile(thisdir, '..', '..');
    addpath(genpath(fullfile(reporoot, 'src')));
end

fprintf('=== Demo: Step 4 — Fitting Infrastructure ===\n');

% generate synthetic problem: intercept + 2 features
n = 200;           % time bins
p = 3;             % intercept + 2 features
rng(777);
X = [ones(n,1), randn(n,1), randn(n,1)];
w_true = [8; 1.5; -0.8];
noise_std = 0.5;
y = X*w_true + noise_std*randn(n,1);

% simple second-difference penalty for illustration
D = sparse(p-2, p);
for i = 1:p-2
    D(i,i) = 1; D(i,i+1) = -2; D(i,i+2) = 1;
end

% lambda grid
lambda_grid = [0, 0.1, 1.0];
fits = cell(numel(lambda_grid),1);

% fit models across lambdas
for i = 1:numel(lambda_grid)
    cfg = struct();
    cfg.model = struct('lambda', lambda_grid(i), 'penalty_operator', D, 'link', 'identity');
    cfg.optim = struct('display','off','max_iter',200,'tol_fun',1e-6,'tol_grad',1e-5);
    fits{i} = fit_glm_map_hg(X, y, cfg);

    fprintf('\n-- lambda = %.3f --\n', lambda_grid(i));
    fprintf('true w:    %8.3f %8.3f %8.3f\n', w_true);
    fprintf('fitted w:  %8.3f %8.3f %8.3f\n', fits{i}.w);
    fprintf('converged: %d | exitflag: %d\n', fits{i}.converged, fits{i}.exitflag);
    nIter = NaN; if isfield(fits{i}.output,'iterations'), nIter = fits{i}.output.iterations; end
    fprintf('iters: %s | grad_norm: %.3e | nll: %.6f\n', mat2str(nIter), fits{i}.grad_norm, fits{i}.nll);
end

% identify best lambda by correlation
corrs = cellfun(@(f) corr(f.mu, y), fits);
[~, bestIdx] = max(corrs);

% plot shrinkage curve and scatter for best lambda
figure('Name','Step 4: Fitting Demo','Color','w');
subplot(1,2,1);
plot(lambda_grid, cellfun(@(f) norm(f.w), fits), '-o','LineWidth',1.5);
xlabel('\lambda'); ylabel('||w||_2'); title('Shrinkage vs \lambda'); grid on;

subplot(1,2,2);
scatter(y, fits{bestIdx}.mu, 15, 'filled'); hold on;
lsline; xlabel('Observed y'); ylabel('Predicted \mu');
title(sprintf('Obs vs Pred (best \lambda=%.2g), r=%.3f', lambda_grid(bestIdx), corrs(bestIdx)));
grid on;

% simple success prints
fprintf('\nBest lambda (by corr): %.3f | r=%.3f\n', lambda_grid(bestIdx), corrs(bestIdx));
fprintf('\n\xE2\x9C\x93 Step 4 Complete: Fitting infrastructure working\n');

