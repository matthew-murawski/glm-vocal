function out = metrics(y, mu)
%METRICS Placeholder computation of evaluation metrics.
%   out = METRICS(y, mu) returns zeroed summary fields.

unused = {y, mu}; %#ok<NASGU>
out = struct('nll', 0, 'pseudoR2', 0, 'deviance_explained', 0);
end
