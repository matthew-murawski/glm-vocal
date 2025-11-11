function kernels = extract_kernels(w, colmap, stim)
%EXTRACT_KERNELS Parse fitted weights and extract kernel structures.
%   kernels = EXTRACT_KERNELS(w, colmap, stim) extracts and reshapes
%   fitted weight vector into interpretable kernel structures for each
%   regressor type. Includes time vectors for plotting.
%
%   Inputs:
%       w      - [p × 1] fitted weight vector
%       colmap - struct from assemble_design_matrix with column mappings
%       stim   - struct with dt field for time vector construction
%
%   Output:
%       kernels - struct with fields for each regressor type:
%                 * intercept: scalar intercept value
%                 * heard_<type>: struct with weights, times, info
%                 * produced_<type>: struct with weights, times, info
%                 * states: struct with convo and spon coefficients
%
%   Each kernel struct contains:
%       * weights: [n_lags × 1] kernel weights
%       * times: [n_lags × 1] time vector (seconds)
%       * info: metadata from colmap (window, mode, basis info)

% validate inputs
if ~isvector(w)
    error('glm:InvalidInput', 'w must be a vector');
end

w = w(:);  % ensure column vector

if ~isstruct(colmap)
    error('glm:InvalidInput', 'colmap must be a struct');
end

if ~isstruct(stim) || ~isfield(stim, 'dt')
    error('glm:InvalidInput', 'stim must be a struct with dt field');
end

dt = stim.dt;

% initialize output struct
kernels = struct();

% section extract intercept
if isfield(colmap, 'intercept')
    intercept_col = colmap.intercept.cols;
    if intercept_col <= length(w)
        kernels.intercept = w(intercept_col);
    end
end

% section extract heard kernels
if isfield(colmap, 'heard_fields') && ~isempty(colmap.heard_fields)
    for i = 1:length(colmap.heard_fields)
        field = colmap.heard_fields{i};

        if isfield(colmap, field)
            cols = colmap.(field).cols;

            % validate column indices
            if any(cols < 1) || any(cols > length(w))
                warning('glm:InvalidCols', ...
                    'Column indices for %s out of bounds. Skipping.', field);
                continue;
            end

            % extract basis weights
            basis_weights = w(cols);

            % extract info
            if isfield(colmap.(field), 'info')
                info = colmap.(field).info;

                % check if basis projection was used
                if isfield(info, 'basis') && ~isempty(info.basis) && ...
                        isfield(info.basis, 'matrix') && ~isempty(info.basis.matrix)
                    % reconstruct full kernel from basis weights
                    basis_matrix = info.basis.matrix;
                    kernel_weights = basis_matrix * basis_weights;

                    % use lag times from info
                    if isfield(info, 'lag_times_s') && ~isempty(info.lag_times_s)
                        kernel_times = info.lag_times_s(:);
                    else
                        % shouldn't happen, but create from lag bins
                        kernel_times = info.lag_bins * dt;
                    end
                else
                    % no basis projection, weights are the kernel
                    kernel_weights = basis_weights;

                    % build time vector
                    if isfield(info, 'lag_times_s') && ~isempty(info.lag_times_s)
                        kernel_times = info.lag_times_s(:);
                    elseif isfield(info, 'lag_bins') && ~isempty(info.lag_bins)
                        kernel_times = info.lag_bins * dt;
                    else
                        % fallback
                        kernel_times = ((0:length(kernel_weights)-1) * dt)';
                    end
                end
            else
                % no info available
                kernel_weights = basis_weights;
                kernel_times = ((0:length(kernel_weights)-1) * dt)';
                info = struct();
            end

            % store in output
            kernels.(field) = struct();
            kernels.(field).weights = kernel_weights(:);
            kernels.(field).times = kernel_times(:);
            kernels.(field).info = info;
        end
    end
end

% section extract produced kernels
if isfield(colmap, 'produced_fields') && ~isempty(colmap.produced_fields)
    for i = 1:length(colmap.produced_fields)
        field = colmap.produced_fields{i};

        if isfield(colmap, field)
            cols = colmap.(field).cols;

            % validate column indices
            if any(cols < 1) || any(cols > length(w))
                warning('glm:InvalidCols', ...
                    'Column indices for %s out of bounds. Skipping.', field);
                continue;
            end

            % extract basis weights
            basis_weights = w(cols);

            % extract info and build time vector (same logic as heard)
            if isfield(colmap.(field), 'info')
                info = colmap.(field).info;

                % check if basis projection was used
                if isfield(info, 'basis') && ~isempty(info.basis) && ...
                        isfield(info.basis, 'matrix') && ~isempty(info.basis.matrix)
                    % reconstruct full kernel from basis weights
                    basis_matrix = info.basis.matrix;
                    kernel_weights = basis_matrix * basis_weights;

                    % use lag times from info
                    if isfield(info, 'lag_times_s') && ~isempty(info.lag_times_s)
                        kernel_times = info.lag_times_s(:);
                    else
                        kernel_times = info.lag_bins * dt;
                    end
                else
                    % no basis projection
                    kernel_weights = basis_weights;

                    if isfield(info, 'lag_times_s') && ~isempty(info.lag_times_s)
                        kernel_times = info.lag_times_s(:);
                    elseif isfield(info, 'lag_bins') && ~isempty(info.lag_bins)
                        kernel_times = info.lag_bins * dt;
                    else
                        kernel_times = ((0:length(kernel_weights)-1) * dt)';
                    end
                end
            else
                kernel_weights = basis_weights;
                kernel_times = ((0:length(kernel_weights)-1) * dt)';
                info = struct();
            end

            % store in output
            kernels.(field) = struct();
            kernels.(field).weights = kernel_weights(:);
            kernels.(field).times = kernel_times(:);
            kernels.(field).info = info;
        end
    end
end

% section extract state coefficients
if isfield(colmap, 'states')
    state_info = colmap.states;

    if isfield(state_info, 'cols')
        cols = state_info.cols;

        % validate column indices
        if all(cols >= 1 & cols <= length(w))
            kernels.states = struct();

            % extract individual state coefficients
            if isfield(state_info, 'names') && length(state_info.names) == length(cols)
                for i = 1:length(state_info.names)
                    state_name = state_info.names{i};
                    kernels.states.(state_name) = w(cols(i));
                end
            else
                % fallback: just store all state weights
                kernels.states.weights = w(cols);
            end
        end
    end
end

end
