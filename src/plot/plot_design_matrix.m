function plot_design_matrix(X, colmap, opts)
% section defaults and validation
% prepare plotting options, coerce the matrix to sparse, and exit early when nothing is selected.
if nargin < 3 || isempty(opts)
    opts = struct();
end
if isempty(X)
    return
end

if issparse(X)
    Xplot = X;
else
    Xplot = sparse(X);
end

nCols = size(Xplot, 2);
if nCols == 0
    return
end

if ~isfield(opts, 'col_range') || isempty(opts.col_range)
    maxCols = min(50, nCols);
    cols = 1:maxCols;
else
    cols = opts.col_range(:)';
end

cols = cols(cols >= 1 & cols <= nCols);
if isempty(cols)
    return
end

outputPath = '';
if isfield(opts, 'output_path') && ~isempty(opts.output_path)
    outputPath = char(opts.output_path);
end

% section figure setup
% generate a grayscale heatmap of the selected columns with optional labels derived from the column map metadata.
fig = figure('Visible', 'off');
try
    data = full(Xplot(:, cols));
    imagesc(data);
    colormap(fig, 'gray');
    colorbar;
    xlabel('regressor columns');
    ylabel('time bins');
    title('design matrix preview');

    labels = makeColumnLabels(cols, colmap);
    set(gca, 'XTick', 1:numel(cols));
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 45);

    if ~isempty(outputPath)
        exportgraphics(fig, outputPath, 'ContentType', 'vector');
    end
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end

function labels = makeColumnLabels(cols, colmap)
% section label helper
% derive human-friendly labels for the selected columns using the column map when available.
labels = cell(size(cols));
for ii = 1:numel(cols)
    col = cols(ii);
    labels{ii} = sprintf('col%d', col);
end

if nargin < 2 || isempty(colmap)
    return
end

fields = fieldnames(colmap);
for ff = 1:numel(fields)
    entries = colmap.(fields{ff});
    if ~isstruct(entries)
        continue
    end
    for ee = 1:numel(entries)
        entry = entries(ee);
        if ~isfield(entry, 'cols')
            continue
        end
        entryCols = entry.cols;
        for ii = 1:numel(cols)
            match = find(entryCols == cols(ii), 1);
            if isempty(match)
                continue
            end

            if isfield(entry, 'names')
                namesField = entry.names;
                if isstring(namesField)
                    namesField = cellstr(namesField);
                elseif ~iscell(namesField)
                    namesField = {namesField};
                end
            else
                namesField = {};
            end

            if numel(namesField) >= match
                labels{ii} = namesField{match};
            else
                labels{ii} = fields{ff};
            end
        end
    end
end
end
