function [wmap, fitinfo] = fit_glm_map(Xd, D, lambda, optCfg)
%FIT_GLM_MAP Placeholder MAP fitter for the GLM.
%   [wmap, fitinfo] = FIT_GLM_MAP(Xd, D, lambda, optCfg) returns empty results.

unused = {Xd, D, lambda, optCfg}; %#ok<NASGU>
wmap = struct('w', zeros(0, 1));
fitinfo = struct();
end
