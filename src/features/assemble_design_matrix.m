function Xd = assemble_design_matrix(streams, states, sps, cfg)
%ASSEMBLE_DESIGN_MATRIX Placeholder assembler for the design matrix.
%   Xd = ASSEMBLE_DESIGN_MATRIX(streams, states, sps, cfg) returns empty scaffolding.

unused = {streams, states, sps, cfg}; %#ok<NASGU>
Xd = struct('X', sparse([], [], [], 0, 0), 'y', zeros(0, 1), 'colmap', struct());
end
