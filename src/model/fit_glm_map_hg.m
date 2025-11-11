function fit = fit_glm_map_hg(X, y, cfg)
%FIT_GLM_MAP_HG Fit Gaussian GLM for high gamma power data.
%   fit = FIT_GLM_MAP_HG(X, y, cfg) fits a Gaussian GLM with regularization
%   using maximum a posteriori (MAP) estimation via gradient-based optimization.
%
%   Inputs:
%       X   - [n × p] design matrix (can be sparse)
%       y   - [n × 1] observed power trace
%       cfg - configuration struct with fields:
%             * model.lambda: regularization strength (default: 0)
%             * model.penalty_operator: D matrix (default: empty)
%             * model.link: 'identity' (default, only option for now)
%             * optim.max_iter: maximum iterations (default: 200)
%             * optim.tol_fun: function tolerance (default: 1e-6)
%             * optim.tol_grad: gradient tolerance (default: 1e-5)
%             * optim.display: 'off'|'iter'|'final' (default: 'off')
%
%   Output:
%       fit - struct with fields:
%             * w: [p × 1] fitted weights
%             * mu: [n × 1] predicted mean power
%             * nll: final negative log-likelihood
%             * exitflag: convergence status from optimizer
%             * output: optimizer output struct
%             * converged: boolean (true if exitflag > 0)
%             * grad_norm: final gradient magnitude

% set defaults for missing config fields
if ~isfield(cfg, 'model')
    cfg.model = struct();
end
if ~isfield(cfg.model, 'lambda')
    cfg.model.lambda = 0;
end
if ~isfield(cfg.model, 'penalty_operator')
    cfg.model.penalty_operator = sparse(0, size(X, 2));
end
if ~isfield(cfg.model, 'link')
    cfg.model.link = 'identity';
end

if ~isfield(cfg, 'optim')
    cfg.optim = struct();
end
if ~isfield(cfg.optim, 'max_iter')
    cfg.optim.max_iter = 200;
end
if ~isfield(cfg.optim, 'tol_fun')
    cfg.optim.tol_fun = 1e-6;
end
if ~isfield(cfg.optim, 'tol_grad')
    cfg.optim.tol_grad = 1e-5;
end
if ~isfield(cfg.optim, 'display')
    cfg.optim.display = 'off';
end

% validate inputs
if ~ismatrix(X)
    error('glm:InvalidInput', 'X must be a matrix.');
end
if ~isvector(y)
    error('glm:InvalidInput', 'y must be a vector.');
end

% extract parameters
lambda = cfg.model.lambda;
D = cfg.model.penalty_operator;
link = cfg.model.link;

% validate link function
if ~strcmp(link, 'identity')
    error('glm:InvalidLink', 'Only identity link is supported in this step. Got: %s', link);
end

% get dimensions
[n, p] = size(X);
if length(y) ~= n
    error('glm:DimensionMismatch', 'X has %d rows but y has %d elements.', n, length(y));
end

% ensure y is column vector
y = y(:);

% initialize weights (intercept-only model)
w0 = zeros(p, 1);

% create objective function handle
objFun = @(w) neglogli_gaussian(w, X, y, lambda, D);

% set up optimization options
try
    % newer MATLAB syntax (R2016a+)
    options = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...
        'Display', cfg.optim.display, ...
        'MaxIterations', cfg.optim.max_iter, ...
        'MaxFunctionEvaluations', cfg.optim.max_iter * 10, ...
        'OptimalityTolerance', cfg.optim.tol_grad, ...
        'StepTolerance', cfg.optim.tol_fun, ...
        'SpecifyObjectiveGradient', true);
catch
    % older MATLAB syntax (pre-R2016a)
    options = optimset(...
        'GradObj', 'on', ...
        'Display', cfg.optim.display, ...
        'MaxIter', cfg.optim.max_iter, ...
        'MaxFunEvals', cfg.optim.max_iter * 10, ...
        'TolFun', cfg.optim.tol_fun, ...
        'TolX', cfg.optim.tol_grad);
end

% run optimization
if exist('fminunc', 'file') == 2
    % use fminunc if available (Optimization Toolbox)
    [w_fit, nll_final, exitflag, output] = fminunc(objFun, w0, options);
else
    % fallback to fminsearch (base MATLAB, slower, no gradient)
    warning('glm:NoOptimToolbox', ...
        'fminunc not found. Using fminsearch (slower, ignores gradient).');
    optset = optimset('Display', cfg.optim.display, ...
        'MaxIter', cfg.optim.max_iter, ...
        'TolFun', cfg.optim.tol_fun);
    [w_fit, nll_final, exitflag, output] = fminsearch(@(w) objFun(w), w0, optset);
    output.gradient = [];
end

% compute predictions
mu_fit = predict_power(w_fit, X, link);

% get final gradient norm
if isfield(output, 'gradient') && ~isempty(output.gradient)
    grad_norm = norm(output.gradient);
else
    % recompute gradient if not available
    [~, grad_final] = objFun(w_fit);
    grad_norm = norm(grad_final);
    output.gradient = grad_final;
end

% build output struct
fit = struct();
fit.w = w_fit;
fit.mu = mu_fit;
fit.nll = nll_final;
fit.exitflag = exitflag;
fit.output = output;
fit.converged = (exitflag > 0);
fit.grad_norm = grad_norm;

% check convergence and issue warnings if needed
check_convergence(fit, cfg);

end

%% helper functions

function check_convergence(fit, cfg) %#ok<INUSD>
%CHECK_CONVERGENCE Issue warnings for convergence problems.

if ~fit.converged
    warning('glm:NoConvergence', ...
        'Optimizer did not converge (exitflag=%d). Results may be unreliable.', ...
        fit.exitflag);

    if isfield(fit.output, 'message')
        warning('glm:OptimizerMessage', 'Optimizer message: %s', fit.output.message);
    end
end

% warn if gradient norm is suspiciously high
if fit.grad_norm > 1e-2
    warning('glm:HighGradient', ...
        'Final gradient norm is high (%.4f). Solution may not be at optimum.', ...
        fit.grad_norm);
end

% warn if NLL is not finite
if ~isfinite(fit.nll)
    warning('glm:InvalidNLL', ...
        'Final NLL is %s. Model fitting failed.', ...
        mat2str(fit.nll));
end

end
