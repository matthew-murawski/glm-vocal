function summary = qc_session_summary(Xd, wmap, rate, cvinfo, kernels, ptest, outdir)
% section setup
% generate summary artifacts (figure, text, json) that capture key qc information for a fitted session.
if nargin < 5 || isempty(kernels)
    kernels = struct();
end
if nargin < 6 || isempty(ptest)
    ptest = struct();
end
if nargin < 7 || isempty(outdir)
    outdir = pwd;
end
if isstring(outdir) || ischar(outdir)
    outdir = char(outdir);
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

summary = struct();
summary.paths = struct();

% section derived stats
% compute basic scalar summaries for downstream reporting.
summary.stats = struct();
summary.stats.n_bins = size(Xd.X, 1);
summary.stats.n_params = numel(wmap.w);
summary.stats.weight_norm = norm(wmap.w);
if isfield(cvinfo, 'lambdas') && isfield(cvinfo, 'mean_nll')
    [minVal, idx] = min(cvinfo.mean_nll);
    summary.stats.best_lambda = cvinfo.lambdas(idx);
    summary.stats.best_lambda_mean_nll = minVal;
else
    summary.stats.best_lambda = NaN;
    summary.stats.best_lambda_mean_nll = NaN;
end
if isfield(rate, 'metrics')
    summary.stats.metrics = rate.metrics;
end

summary.stats.event_counts = normalize_event_counts(rate);
if isfield(kernels, 'states') && isstruct(kernels.states)
    summary.stats.state_coeffs = struct();
    if isfield(kernels.states, 'convo')
        summary.stats.state_coeffs.convo = kernels.states.convo;
    end
    if isfield(kernels.states, 'spon')
        summary.stats.state_coeffs.spon = kernels.states.spon;
    end
end

summary.stats.permutation = summarize_permutation_tests(ptest);

% section statistical test
% perform a wald test for the difference between conversational and spontaneous state coefficients.
if isfield(wmap, 'hessian') && isfield(Xd, 'colmap') && isfield(Xd.colmap, 'states')
    colmap = Xd.colmap.states;
    if isfield(colmap, 'convo') && isfield(colmap, 'spon')
        try
            cov_matrix = inv(wmap.hessian);
            idx_convo = colmap.convo;
            idx_spon = colmap.spon;

            w_convo = wmap.w(idx_convo);
            w_spon = wmap.w(idx_spon);
            var_convo = cov_matrix(idx_convo, idx_convo);
            var_spon = cov_matrix(idx_spon, idx_spon);
            cov_convo_spon = cov_matrix(idx_convo, idx_spon);

            diff = w_convo - w_spon;
            se_diff = sqrt(var_convo + var_spon - 2 * cov_convo_spon);
            z_score = diff / se_diff;
            p_value = 2 * (1 - normcdf(abs(z_score)));

            summary.stats.state_coeffs.difference = diff;
            summary.stats.state_coeffs.z_score = z_score;
            summary.stats.state_coeffs.p_value = p_value;
        catch err
            warning('qc_session_summary:WaldTestFailed', 'could not compute wald test: %s', err.message);
        end
    end
end


% section summary figure
% compose a lightweight summary figure with rate traces, kernel weights, and cv curve where available.
fig = figure('Visible', 'off');
try
    tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');

    % subplot 1: rate vs spikes if available
    nexttile;
    if isfield(rate, 'stim') && isfield(rate.stim, 't') && isfield(rate, 'y') && isfield(rate, 'mu')
        t = rate.stim.t(:);
        y = rate.y(:);
        mu = rate.mu(:);
        plot(t, y, 'Color', [0.3 0.3 0.3], 'DisplayName', 'spikes');
        hold on
        plot(t, mu, 'r', 'LineWidth', 1.2, 'DisplayName', 'rate');
        hold off
        xlabel('time (s)'); ylabel('count');
        title('rate vs spikes');
        legend('Location', 'best');
    else
        text(0.5, 0.5, 'rate data unavailable', 'HorizontalAlignment', 'center');
        axis off
    end

    % subplot 2: weight stem plot
    nexttile;
    stem(wmap.w, 'Marker', 'none');
    xlabel('parameter index'); ylabel('weight');
    title('weights'); grid on;

    % subplot 3: cv curve if provided
    nexttile;
    if isfield(cvinfo, 'lambdas') && isfield(cvinfo, 'mean_nll')
        semilogx(cvinfo.lambdas, cvinfo.mean_nll, '-o', 'LineWidth', 1.2);
        xlabel('\lambda'); ylabel('mean nll');
        title('cv curve'); grid on;
    else
        text(0.5, 0.5, 'cv info unavailable', 'HorizontalAlignment', 'center');
        axis off
    end

    % subplot 4: textual summary
    nexttile;
    axis off;
    textLines = {
        sprintf('bins: %d', summary.stats.n_bins), ...
        sprintf('params: %d', summary.stats.n_params), ...
        sprintf('||w||: %.3f', summary.stats.weight_norm)
    };
    textLines{end+1} = ''; %#ok<AGROW>
    textLines{end+1} = 'Event Counts:'; %#ok<AGROW>
    textLines{end+1} = format_count_line('heard', summary.stats.event_counts.heard); %#ok<AGROW>
    if ~isnan(summary.stats.event_counts.produced_any)
        textLines{end+1} = format_count_line('produced any', summary.stats.event_counts.produced_any); %#ok<AGROW>
    end
    producedFields = summary.stats.event_counts.produced_fields;
    if isempty(producedFields)
        textLines{end+1} = '  produced: n/a'; %#ok<AGROW>
    else
        for ii = 1:numel(producedFields)
            fname = producedFields{ii};
            label = format_produced_label(fname);
            value = summary.stats.event_counts.produced.(fname);
            textLines{end+1} = format_count_line(label, value); %#ok<AGROW>
        end
    end
    if ~isnan(summary.stats.best_lambda)
        textLines{end+1} = sprintf('best lambda: %.3g', summary.stats.best_lambda); %#ok<AGROW>
        textLines{end+1} = sprintf('mean nll: %.3f', summary.stats.best_lambda_mean_nll); %#ok<AGROW>
    end
    if isfield(summary.stats, 'state_coeffs')
        textLines{end+1} = ''; %#ok<AGROW>
        textLines{end+1} = 'State Coeffs:'; %#ok<AGROW>
        if isfield(summary.stats.state_coeffs, 'convo')
            textLines{end+1} = sprintf('  convo: %.4f', summary.stats.state_coeffs.convo); %#ok<AGROW>
        end
        if isfield(summary.stats.state_coeffs, 'spon')
            textLines{end+1} = sprintf('  spon: %.4f', summary.stats.state_coeffs.spon); %#ok<AGROW>
        end
        if isfield(summary.stats.state_coeffs, 'difference')
            textLines{end+1} = sprintf('  diff (convo-spon): %.4f', summary.stats.state_coeffs.difference); %#ok<AGROW>
            textLines{end+1} = sprintf('  z=%.2f, p=%.3f', summary.stats.state_coeffs.z_score, summary.stats.state_coeffs.p_value); %#ok<AGROW>
        end
    end
    if summary.stats.permutation.n_kernels > 0
        textLines{end+1} = ''; %#ok<AGROW>
        textLines{end+1} = sprintf('Permutation (p<%.2f):', summary.stats.permutation.threshold); %#ok<AGROW>
        textLines{end+1} = sprintf('  significant: %d / %d', summary.stats.permutation.n_significant, summary.stats.permutation.n_kernels); %#ok<AGROW>
        permNames = fieldnames(summary.stats.permutation.kernels);
        for ii = 1:min(4, numel(permNames))
            fname = permNames{ii};
            entry = summary.stats.permutation.kernels.(fname);
            textLines{end+1} = sprintf('  %s: p=%.3f', fname, entry.p_value); %#ok<AGROW>
        end
        if numel(permNames) > 4
            textLines{end+1} = sprintf('  (+%d more)', numel(permNames) - 4); %#ok<AGROW>
        end
    end
    text(0, 1, strjoin(textLines, newline), 'VerticalAlignment', 'top');

    figPath = fullfile(outdir, 'qc_summary.png');
    saveas(fig, figPath);
    summary.paths.figure = figPath;
catch err
    close(fig);
    rethrow(err);
end
close(fig);

% section text summary
% write a human-readable text file summarising key stats.
textPath = fullfile(outdir, 'qc_summary.txt');
fid = fopen(textPath, 'w');
if fid ~= -1
    nl = newline;
    fprintf(fid, 'QC SUMMARY%s', nl);
    fprintf(fid, 'bins: %d%s', summary.stats.n_bins, nl);
    fprintf(fid, 'params: %d%s', summary.stats.n_params, nl);
    fprintf(fid, '||w||: %.4f%s', summary.stats.weight_norm, nl);
    if ~isnan(summary.stats.best_lambda)
        fprintf(fid, 'best lambda: %.6g%s', summary.stats.best_lambda, nl);
        fprintf(fid, 'mean nll (best): %.4f%s', summary.stats.best_lambda_mean_nll, nl);
    end
    fprintf(fid, '%sEvent Counts:%s', nl, nl);
    fprintf(fid, '%s%s', format_count_text('heard', summary.stats.event_counts.heard), nl);
    if ~isnan(summary.stats.event_counts.produced_any)
        fprintf(fid, '%s%s', format_count_text('produced any', summary.stats.event_counts.produced_any), nl);
    end
    producedFields = summary.stats.event_counts.produced_fields;
    if isempty(producedFields)
        fprintf(fid, '  produced: n/a%s', nl);
    else
        for ii = 1:numel(producedFields)
            fname = producedFields{ii};
            label = format_produced_label(fname);
            fprintf(fid, '%s%s', format_count_text(label, summary.stats.event_counts.produced.(fname)), nl);
        end
    end
    if isfield(summary.stats, 'state_coeffs')
        fprintf(fid, '%sState Coefficients:%s', nl, nl);
        if isfield(summary.stats.state_coeffs, 'convo')
            fprintf(fid, '  convo: %.4f%s', summary.stats.state_coeffs.convo, nl);
        end
        if isfield(summary.stats.state_coeffs, 'spon')
            fprintf(fid, '  spon: %.4f%s', summary.stats.state_coeffs.spon, nl);
        end
        if isfield(summary.stats.state_coeffs, 'difference')
            fprintf(fid, '  diff (convo-spon): %.4f%s', summary.stats.state_coeffs.difference, nl);
            fprintf(fid, '  z-score: %.2f%s', summary.stats.state_coeffs.z_score, nl);
            fprintf(fid, '  p-value: %.3f%s', summary.stats.state_coeffs.p_value, nl);
        end
    end
    if summary.stats.permutation.n_kernels > 0
        fprintf(fid, '%sPermutation Tests (threshold %.2f):%s', nl, summary.stats.permutation.threshold, nl);
        fprintf(fid, '  significant: %d / %d%s', summary.stats.permutation.n_significant, summary.stats.permutation.n_kernels, nl);
        permNames = fieldnames(summary.stats.permutation.kernels);
        for ii = 1:numel(permNames)
            fname = permNames{ii};
            entry = summary.stats.permutation.kernels.(fname);
            fprintf(fid, '  %s: p=%.3f%s', fname, entry.p_value, nl);
        end
    end
    fclose(fid);
else
    warning('qc_session_summary:FileWriteFailed', 'unable to write text summary to %s', textPath);
end
summary.paths.text = textPath;

% section json summary
% encode the summary struct to json for downstream tooling.
jsonPath = fullfile(outdir, 'qc_summary.json');
try
    jsonStr = jsonencode(summary);
    fid = fopen(jsonPath, 'w');
    if fid ~= -1
        fwrite(fid, jsonStr, 'char');
        fclose(fid);
    else
        warning('qc_session_summary:FileWriteFailed', 'unable to write json summary to %s', jsonPath);
    end
catch err
    warning('qc_session_summary:JsonEncodeFailed', 'json encoding failed: %s', err.message);
end
summary.paths.json = jsonPath;
end

function permSummary = summarize_permutation_tests(ptest)
% section permutation helper
% convert raw permutation test outputs into a compact summary that highlights significant kernels.
permSummary = struct('threshold', 0.05, 'n_kernels', 0, 'n_significant', 0, 'kernels', struct());

if ~isstruct(ptest) || isempty(fieldnames(ptest))
    return
end

fields = fieldnames(ptest);
for ii = 1:numel(fields)
    fname = fields{ii};
    entry = ptest.(fname);
    if ~isstruct(entry) || ~isfield(entry, 'p_value')
        continue
    end

    pVal = entry.p_value;
    if ~isscalar(pVal) || ~isfinite(pVal)
        continue
    end

    permSummary.n_kernels = permSummary.n_kernels + 1;
    isSig = pVal < permSummary.threshold;

    kernelSummary = struct('p_value', pVal, 'significant', isSig);
    if isfield(entry, 'ci_lower') && isfield(entry, 'ci_upper')
        kernelSummary.ci_lower = entry.ci_lower(:)';
        kernelSummary.ci_upper = entry.ci_upper(:)';
    end
    permSummary.kernels.(fname) = kernelSummary;

    if isSig
        permSummary.n_significant = permSummary.n_significant + 1;
    end
end
end

function line = format_count_line(label, value)
% section figure text helper
% format event counts for inclusion in the summary figure.
if isnan(value)
    line = sprintf('  %s: n/a', label);
else
    line = sprintf('  %s: %d', label, round(value));
end
end

function textLine = format_count_text(label, value)
% section text helper
% format event counts for the plain-text summary file.
if isnan(value)
    textLine = sprintf('  %s: n/a', label);
else
    textLine = sprintf('  %s: %d', label, round(value));
end
end

function counts = normalize_event_counts(rate)
counts = struct('heard', NaN, 'produced_fields', {{}}, 'produced', struct(), 'produced_any', NaN);
if ~isstruct(rate) || ~isfield(rate, 'event_counts') || ~isstruct(rate.event_counts)
    return
end

raw = rate.event_counts;
if isfield(raw, 'heard') && isscalar(raw.heard)
    counts.heard = double(raw.heard);
end
if isfield(raw, 'produced_any') && isscalar(raw.produced_any)
    counts.produced_any = double(raw.produced_any);
end

if isfield(raw, 'produced_fields') && ~isempty(raw.produced_fields)
    fields = cellstr(raw.produced_fields(:));
    counts.produced_fields = fields;
    for ii = 1:numel(fields)
        fname = fields{ii};
        if isfield(raw, 'produced') && isfield(raw.produced, fname)
            counts.produced.(fname) = double(raw.produced.(fname));
        else
            counts.produced.(fname) = NaN;
        end
    end
else
    legacy = {'produced_spontaneous', 'produced_after_heard', 'produced_after_produced'};
    available = legacy(isfield(raw, legacy));
    counts.produced_fields = available;
    for ii = 1:numel(available)
        fname = available{ii};
        counts.produced.(fname) = double(raw.(fname));
    end
    if isnan(counts.produced_any) && ~isempty(available)
        vals = cellfun(@(f) counts.produced.(f), available);
        if all(isfinite(vals))
            counts.produced_any = sum(vals);
        end
    end
end
end

function label = format_produced_label(fieldName)
label = strrep(fieldName, 'produced_', 'produced ');
label = strrep(label, '_', ' ');
end
