function summary = qc_session_summary(Xd, wmap, rate, cvinfo, outdir)
% section setup
% generate summary artifacts (figure, text, json) that capture key qc information for a fitted session.
if nargin < 5 || isempty(outdir)
    outdir = pwd;
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
    if ~isnan(summary.stats.best_lambda)
        textLines{end+1} = sprintf('best lambda: %.3g', summary.stats.best_lambda); %#ok<AGROW>
        textLines{end+1} = sprintf('mean nll: %.3f', summary.stats.best_lambda_mean_nll); %#ok<AGROW>
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
