function plot_cv_curve(cvinfo, outdir)
% section setup
% visualise cross-validation mean nll across lambdas and save figure for qc.
if nargin < 2 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

lambdas = cvinfo.lambdas(:);
meanNll = cvinfo.mean_nll(:);

fig = figure('Visible', 'off');
try
    semilogx(lambdas, meanNll, '-o', 'LineWidth', 1.2);
    xlabel('\lambda');
    ylabel('mean held-out nll per bin');
    title('cv curve');
    grid on;

    if isfield(cvinfo, 'fold_nll')
        hold on
        for k = 1:size(cvinfo.fold_nll, 2)
            semilogx(lambdas, cvinfo.fold_nll(:, k), ':', 'Color', [0.6 0.6 0.6]);
        end
        hold off
    end

    filepath = fullfile(outdir, 'cv_curve.png');
    saveas(fig, filepath);
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end
