function plot_rate_vs_spikes(stim, y, mu, outdir)
% section setup
% plot observed spike counts against predicted rate trajectory and save figure for qc.
if nargin < 4 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

t = stim.t(:);
if numel(t) ~= numel(y) || numel(y) ~= numel(mu)
    error('plot_rate_vs_spikes:SizeMismatch', 'stim, y, and mu must align.');
end

y = y(:);
mu = mu(:);

fig = figure('Visible', 'off');
try
    plot(t, y, 'Color', [0.2 0.2 0.2], 'DisplayName', 'spike count');
    hold on
    plot(t, mu, 'r', 'LineWidth', 1.2, 'DisplayName', 'predicted rate');
    hold off
    xlabel('time (s)');
    ylabel('spikes / rate');
    title('rate vs spikes');
    legend('Location', 'best');
    grid on;

    filepath = fullfile(outdir, 'rate_vs_spikes.png');
    saveas(fig, filepath);
catch err
    close(fig);
    rethrow(err);
end
close(fig);
end
