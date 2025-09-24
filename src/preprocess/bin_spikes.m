function sps = bin_spikes(spike_times, stim)
%BIN_SPIKES Bin spike times onto the provided stimulus grid.
%   sps = BIN_SPIKES(spike_times, stim) returns spike counts per time bin
%   using the time vector and step size stored in stim.

% short-circuit when no spikes are provided
if isempty(spike_times)
    sps = zeros(numel(stim.t), 1);
    return
end

% ensure we are working with a column vector of double times
spike_times = spike_times(:);

% assemble histogram edges to cover the entire stimulus grid
edges = [stim.t; stim.t(end) + stim.dt];

% use histcounts to tally spikes per bin and return a column vector
counts = histcounts(spike_times, edges);
sps = counts(:);
end
