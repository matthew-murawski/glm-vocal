function [nll, grad] = neglogli_poiss(w, X, y)
%NEGLOGLI_POISS Placeholder Poisson negative log-likelihood.
%   [nll, grad] = NEGLOGLI_POISS(w, X, y) returns zero cost and gradient.

unused = {X, y}; %#ok<NASGU>
nll = 0;
grad = zeros(size(w));
end
