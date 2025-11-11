function [eta_trace, eta_time] = calculate_eta_from_continuous(continuous_power, continuous_time, event_onsets, eta_window, target_dt)
% calculate_eta_from_continuous compute an event-triggered average trace.
% given a continuous power vector and matching time vector, this function
% resamples snippets around each event onset onto a common time base and
% averages across events, ignoring nan samples that fall outside the data.

%% section: normalize and validate inputs
% ensure column vectors and basic shapes match; let interp1 handle edge samples.
continuous_power = continuous_power(:);
continuous_time = continuous_time(:);
if numel(continuous_power) ~= numel(continuous_time)
    error('calculate_eta_from_continuous:InvalidInput', 'continuous_power and continuous_time must be the same length');
end

% coerce event onsets into a column vector of doubles
if isempty(event_onsets)
    event_onsets = [];
else
    event_onsets = double(event_onsets(:));
end

%% section: build target time vector for the eta
% we construct a fixed grid relative to t=0 defined by eta_window and target_dt.
eta_time = (eta_window(1):target_dt:eta_window(2))';
nt = numel(eta_time);

% early exit for no events: return nan trace aligned to eta_time
if isempty(event_onsets)
    eta_trace = nan(nt, 1);
    return
end

%% section: resample snippets around each onset
% we use linear interpolation onto the target grid; out-of-range samples remain nan.
nEvents = numel(event_onsets);
snippet_matrix = nan(nEvents, nt);

for ii = 1:nEvents
    t_onset = event_onsets(ii);
    target_times_for_snippet = t_onset + eta_time;
    snippet = interp1(continuous_time, continuous_power, target_times_for_snippet, 'linear');
    snippet_matrix(ii, :) = snippet(:).';
end

%% section: average across events, ignoring nans
% compute the column-wise mean to obtain the final eta trace.
eta_trace = mean(snippet_matrix, 1, 'omitnan');
eta_trace = eta_trace(:);

end
