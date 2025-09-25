function [wmap, fitinfo] = fit_glm_map(Xd, D, lambda, optCfg)
% section inputs and defaults
% fit a poisson glm with quadratic penalty using fminunc when available, otherwise fall back to fminunc-like optimisation.
defaults = struct('tol_fun', 1e-6, 'tol_grad', 1e-5, 'max_iter', 200, 'display', 'off');
if nargin < 4 || isempty(optCfg)
    optCfg = struct();
end
optCfg = mergeStructs(defaults, optCfg);

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
    error('fit_glm_map:InvalidLambda', 'lambda must be a non-negative scalar.');
end

L2op = @(w) l2PenaltyOp(w, D, Dt, lambda);
objFun = @(w) neglogli_poiss(w, X, y, L2op);

if exist('fminunc', 'file') == 2
    options = buildFminuncOptions(optCfg);
    problem = struct();
    problem.objective = @(w) objectiveWrapper(objFun, w);
    problem.x0 = w0;
    problem.solver = 'fminunc';
    problem.options = options;
    [bestW, fval, exitflag, output] = fminunc(problem);
else
    options = struct('Display', optCfg.display, ...
        'MaxIter', optCfg.max_iter, ...
        'TolFun', optCfg.tol_fun);
    [bestW, fval, exitflag, output] = fminsearchFallback(objFun, w0, options);
end

if exitflag <= 0
    warning('fit_glm_map:NoConvergence', 'optimizer did not converge: %s', output.message);
end

wmap = struct();
wmap.w = bestW;

fitinfo = struct();
fitinfo.nll = fval;
fitinfo.exitflag = exitflag;
fitinfo.iterations = output.iterations;
fitinfo.funcCount = output.funcCount;
fitinfo.gradient = output.gradient;
end

function val = objectiveWrapper(objFun, w)
% section objective wrapper
% ensure the function handle returns scalar values compatible with fminunc.
w = w(:);
[val, ~] = objFun(w);
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
        'StepTolerance', optCfg.tol_fun);
catch
    options = struct('Display', optCfg.display, ...
        'MaxIter', optCfg.max_iter, ...
        'MaxFunEvals', optCfg.max_iter * 10, ...
        'TolFun', optCfg.tol_fun, ...
        'TolX', optCfg.tol_grad, ...
        'LargeScale', 'off', ...
        'GradObj', 'off');
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
output = struct('iterations', outputStruct.iterations, ...
    'funcCount', outputStruct.funcCount, ...
    'gradient', grad, ...
    'message', 'Fallback fminsearch');
end
