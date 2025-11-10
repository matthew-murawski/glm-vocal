function [Xd_std, scaling_info] = standardize_design_matrix(Xd, exclude_cols)
%STANDARDIZE_DESIGN_MATRIX Z-score normalize design matrix columns.
%   [Xd_std, scaling_info] = STANDARDIZE_DESIGN_MATRIX(Xd, exclude_cols)
%   standardizes each column of the design matrix to have mean 0 and std 1.
%
%   Inputs:
%     Xd          - struct with fields .X (design matrix), .y (response), .colmap
%     exclude_cols - (optional) cell array of predictor names to exclude from
%                    standardization (e.g., {'intercept'})
%
%   Returns:
%     Xd_std       - struct with standardized .X, original .y, and .colmap
%     scaling_info - struct with .means and .stds for each column
%
%   Note: Columns with zero variance are left unchanged to avoid division by zero.

% section input validation
if nargin < 2 || isempty(exclude_cols)
    exclude_cols = {'intercept'};
elseif ischar(exclude_cols)
    exclude_cols = {exclude_cols};
end

X = Xd.X;
y = Xd.y;
colmap = Xd.colmap;

[nRows, nCols] = size(X);

% section identify columns to standardize
% determine which columns should be excluded from standardization
exclude_mask = false(1, nCols);

% check for intercept column
if isfield(colmap, 'intercept')
    intercept_col = colmap.intercept.cols;
    if intercept_col >= 1 && intercept_col <= nCols
        exclude_mask(intercept_col) = true;
    end
end

% check for other excluded predictors
for ii = 1:numel(exclude_cols)
    field_name = exclude_cols{ii};
    if isfield(colmap, field_name) && ~strcmp(field_name, 'intercept')
        if isfield(colmap.(field_name), 'cols')
            cols = colmap.(field_name).cols;
            valid_cols = cols >= 1 & cols <= nCols;
            exclude_mask(cols(valid_cols)) = true;
        end
    end
end

% section compute statistics
% compute mean and std for each column
means = zeros(1, nCols);
stds = ones(1, nCols);

for col = 1:nCols
    if exclude_mask(col)
        % keep original values for excluded columns
        means(col) = 0;
        stds(col) = 1;
    else
        % compute statistics (handle sparse efficiently)
        if issparse(X)
            col_data = full(X(:, col));
        else
            col_data = X(:, col);
        end

        means(col) = mean(col_data);
        stds(col) = std(col_data);

        % avoid division by zero for constant columns
        if stds(col) < 1e-10
            stds(col) = 1;
            warning('standardize_design_matrix:ZeroVariance', ...
                    'Column %d has near-zero variance (std=%.2e), leaving unstandardized.', ...
                    col, stds(col));
        end
    end
end

% section standardization
% apply z-score normalization: (x - mean) / std
X_std = X;
for col = 1:nCols
    if ~exclude_mask(col) && stds(col) > 0
        if issparse(X)
            % for sparse matrices, center and scale efficiently
            X_std(:, col) = (X(:, col) - means(col)) / stds(col);
        else
            X_std(:, col) = (X(:, col) - means(col)) / stds(col);
        end
    end
end

% section output assembly
Xd_std = struct();
Xd_std.X = X_std;
Xd_std.y = y;
Xd_std.colmap = colmap;
if isfield(Xd, 'response_type')
    Xd_std.response_type = Xd.response_type;
end

scaling_info = struct();
scaling_info.means = means;
scaling_info.stds = stds;
scaling_info.excluded_cols = find(exclude_mask);
scaling_info.n_standardized = sum(~exclude_mask);

end
