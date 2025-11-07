function [wmap, fitinfo] = fit_glm_gauss(Xd, D, lambda, optCfg, link)
%FIT_GLM_GAUSS Fit a Gaussian GLM with L2 smoothness penalty.
%   [wmap, fitinfo] = FIT_GLM_GAUSS(Xd, D, lambda, optCfg, link)
%   fits a Gaussian GLM to continuous response data (LFP) with optional
%   L2 penalty for smoothness.
%
%   Inputs:
%     Xd     - struct with fields .X (design matrix) and .y (response vector)
%     D      - penalty matrix (sparse, usually second-difference operator)
%     lambda - penalty strength (scalar, non-negative)
%     optCfg - (optional) optimization config struct with fields:
%              .tol_fun, .tol_grad, .max_iter, .display
%     link   - (optional) link function: 'identity' (default) or 'log'
%
%   Returns:
%     wmap    - struct with fields:
%               .w (fitted weights)
%               .hessian (Hessian of objective at solution)
%     fitinfo - struct with optimization diagnostics

% section inputs and defaults
% fit a gaussian glm with quadratic penalty using fminunc when available,
% otherwise fall back to fminsearch or closed-form solution for identity link.
defaults = struct('tol_fun', 1e-6, 'tol_grad', 1e-5, 'max_iter', 200, 'display', 'off');
if nargin < 4 || isempty(optCfg)
    optCfg = struct();
end
optCfg = mergeStructs(defaults, optCfg);

if nargin < 5 || isempty(link)
    link = 'identity';
end

X = Xd.X;
y = Xd.y;

nParams = size(X, 2);
w0 = zeros(nParams, 1);

if isempty(D)
    D = sparse(0, nParams);
end
if issparse(D)
    Dt = D';
else
    Dt = D';
end

lambda = double(lambda);
if ~isscalar(lambda) || lambda < 0
    error('fit_glm_gauss:InvalidLambda', 'lambda must be a non-negative scalar.');
end

% section optimization strategy
% for identity link, we can use closed-form solution or iterative methods
% for log link, we must use iterative optimization
if strcmpi(link, 'identity') && issparse(X) && lambda > 0
    % use closed-form solution for identity link: w = (X'*X + lambda*D'*D) \ (X'*y)
    % this is typically faster and more accurate than iterative methods
    try
        XtX = X' * X;
        ridge = lambda * (Dt * D);
        A = XtX + ridge;
        b = X' * y;

        % try sparse Cholesky if A is sparse
        if issparse(A)
            bestW = A \ b;
        else
            bestW = A \ b;
        end

        % compute objective value
        residuals = y - X * bestW;
        nll = 0.5 * sum(residuals.^2);
        Dw = D * bestW;
        penalty = lambda * (Dw' * Dw);
        fval = nll + penalty;

        % gradient should be near zero
        grad = X' * (X * bestW - y) + 2 * lambda * (Dt * (D * bestW));

        % use closed-form Hessian
        hessian = XtX + 2 * lambda * (Dt * D);

        wmap = struct();
        wmap.w = bestW;
        wmap.hessian = hessian;

        fitinfo = struct();
        fitinfo.nll = fval;
        fitinfo.exitflag = 1;  % success
        fitinfo.iterations = 0;  % closed-form
        fitinfo.funcCount = 1;
        fitinfo.gradient = grad;
        fitinfo.method = 'closed_form';

        return
    catch ME
        warning('fit_glm_gauss:ClosedFormFailed', 'Closed-form solution failed: %s. Falling back to iterative optimization.', ME.message);
    end
end

% section iterative optimization
% use fminunc or fminsearch for log link or if closed-form fails
L2op = @(w) l2PenaltyOp(w, D, Dt, lambda);
objFun = @(w) neglogli_gauss(w, X, y, L2op, link);

useFminunc = true;
if isstruct(optCfg) && isfield(optCfg, 'use_fminunc') && ~optCfg.use_fminunc
    useFminunc = false;
end

if useFminunc && exist('fminunc', 'file') == 2
    options = buildFminuncOptions(optCfg);
    problem = struct();
    problem.objective = @(w) objectiveWrapper(objFun, w);
    problem.x0 = w0;
    problem.solver = 'fminunc';
    problem.options = options;
    try
        [bestW, fval, exitflag, output, gradOut] = fminunc(problem);
    catch
        [bestW, fval, exitflag, output] = fminunc(problem);
        gradOut = getStructField(output, 'gradient', zeros(size(bestW)));
    end
    if ~isfield(output, 'gradient') || isempty(output.gradient)
        output.gradient = gradOut;
    end
else
    options = struct('Display', optCfg.display, ...
        'MaxIter', optCfg.max_iter, ...
        'TolFun', optCfg.tol_fun);
    [bestW, fval, exitflag, output] = fminsearchFallback(objFun, w0, options);
end

if exitflag <= 0 && isfield(output, 'message') && useFminunc
    warning('fit_glm_gauss:NoConvergence', 'optimizer did not converge: %s', output.message);
end

% section hessian computation
% compute the hessian of the objective function at the solution
switch lower(link)
    case 'identity'
        % H = X'*X + 2*lambda*D'*D
        hessian = X' * X + 2 * lambda * (Dt * D);

    case 'log'
        % H = X' * diag(mu.^2) * X + 2*lambda*D'*D
        % where mu = exp(X*w) and Hessian accounts for log-link nonlinearity
        Xw = X * bestW;
        Xw = max(min(Xw, 50), -50);  % clip
        mu = exp(Xw);
        mu(mu > 1e8) = 1e8;  % clip for stability
        mu_sq = mu.^2;
        W = spdiags(mu_sq, 0, size(X, 1), size(X, 1));
        H_nll = X' * W * X;
        H_ridge = 2 * lambda * (Dt * D);
        hessian = H_nll + H_ridge;
end

wmap = struct();
wmap.w = bestW;
wmap.hessian = hessian;

fitinfo = struct();
fitinfo.nll = fval;
fitinfo.exitflag = exitflag;
fitinfo.iterations = getStructField(output, 'iterations', NaN);
fitinfo.funcCount = getStructField(output, 'funcCount', NaN);
fitinfo.gradient = getStructField(output, 'gradient', zeros(size(bestW)));
fitinfo.method = 'iterative';
end

function [val, grad] = objectiveWrapper(objFun, w)
% section objective wrapper
% ensure the function handle returns scalar values and gradients compatible with fminunc.
w = w(:);
[val, grad] = objFun(w);
val = double(val);
grad = double(grad(:));
end

function [penalty, grad] = l2PenaltyOp(w, D, Dt, lambda)
% section penalty operator
% compute lambda * ||D w||^2 and its gradient contribution.
if lambda == 0 || isempty(D)
    penalty = 0;
    grad = zeros(size(w));
    return
end

Dw = D * w;
penalty = lambda * (Dw' * Dw);
if issparse(Dt)
    grad = 2 * lambda * full(Dt * Dw);
else
    grad = 2 * lambda * (Dt * Dw);
end
end

function s = mergeStructs(defaults, overrides)
% section struct merge
% merge override fields into defaults while keeping unspecified defaults in place.
s = defaults;
if ~isstruct(overrides)
    return
end
fields = fieldnames(overrides);
for ii = 1:numel(fields)
    s.(fields{ii}) = overrides.(fields{ii});
end
end

function options = buildFminuncOptions(optCfg)
% section options builder
% create optimisation options compatible with available toolboxes.
try
    options = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...
        'Display', optCfg.display, ...
        'MaxIterations', optCfg.max_iter, ...
        'MaxFunctionEvaluations', optCfg.max_iter * 10, ...
        'OptimalityTolerance', optCfg.tol_grad, ...
        'StepTolerance', optCfg.tol_fun, ...
        'SpecifyObjectiveGradient', true);
catch
    options = struct('Display', optCfg.display, ...
        'MaxIter', optCfg.max_iter, ...
        'MaxFunEvals', optCfg.max_iter * 10, ...
        'TolFun', optCfg.tol_fun, ...
        'TolX', optCfg.tol_grad, ...
        'LargeScale', 'off', ...
        'GradObj', 'on');
end
end

function [w, fval, exitflag, output] = fminsearchFallback(objFun, w0, options)
% section fallback optimiser
% use fminsearch (available in base MATLAB) when fminunc is unavailable.
optset = optimset('Display', options.Display, 'MaxIter', options.MaxIter, 'TolFun', options.TolFun, 'TolX', options.TolFun);

wrapper = @(w) objectiveWrapper(objFun, w);
[wFlat, fval, exitflag, outputStruct] = fminsearch(wrapper, w0(:)', optset);
w = wFlat(:);

% fminsearch does not provide gradient, so approximate via objFun
[~, grad] = objFun(w);
output = struct('iterations', getStructField(outputStruct, 'iterations', NaN), ...
    'funcCount', getStructField(outputStruct, 'funcCount', NaN), ...
    'gradient', grad, ...
    'message', 'Fallback fminsearch');
end

function val = getStructField(s, fieldName, defaultVal)
% section field helper
% retrieve a struct field if it exists; otherwise return a default value.
if isstruct(s) && isfield(s, fieldName)
    val = s.(fieldName);
else
    val = defaultVal;
end
end
