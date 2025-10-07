function h_fig = plot_event_prediction_preview(result, options)
% PLOT_EVENT_PREDICTION_PREVIEW Visualizes model performance around a specific event.
%
%   h_fig = plot_event_prediction_preview(result, 'Name', Value, ...)
%
%   Generates a plot showing the binned spikes, model-predicted firing rate,
%   and event markers for a specified time window centered on a particular
%   behavioral event. This is useful for quick quality control of a GLM fit.
%
%   Inputs:
%       result    (struct) - The output struct from run_fit_single_neuron. Must
%                            contain fields: events, rate.mu, sps, and stim.
%
%   Optional Name-Value Pair Inputs:
%       'EventType'   (string) - The type of event to center on (e.g., 'produced').
%                                Default: 'produced'.
%       'EventIndex'  (integer)- Which event of the specified type to use (e.g.,
%                                the 1st, 2nd...). Default: 1.
%       'Window'      (1x2 double) - Time window [before, after] in seconds.
%                                Default: [5, 5].
%       'Title'       (string) - Custom title for the plot. If empty, a default
%                                title is generated. Default: "".
%
%   Outputs:
%       h_fig       (handle) - The handle to the generated figure. Returns empty
%                              if the plot could not be created.

arguments
    result (1,1) struct
    options.EventType   (1,1) string = "produced"
    options.EventIndex  (1,1) double {mustBeInteger, mustBePositive} = 5
    options.Window      (1,2) double {mustBeNumeric} = [10, 10]
    options.Title       (1,1) string = ""
end

h_fig = []; % Initialize output to empty

% --- 1. Find the target event ---
event_mask = arrayfun(@(evt) isfield(evt, 'kind') && strcmpi(evt.kind, options.EventType), result.events);
target_events = result.events(event_mask);

if isempty(target_events)
    warning('No events of type "%s" found. Cannot generate plot.', options.EventType);
    return;
end

if options.EventIndex > numel(target_events)
    warning('Requested event index %d, but only %d events of type "%s" found.', ...
        options.EventIndex, numel(target_events), options.EventType);
    return;
end

target_event = target_events(options.EventIndex);
event_time_sec = target_event.t_on;

if isempty(event_time_sec) || ~isfinite(event_time_sec)
    warning('Invalid timestamp for the selected event. Cannot generate plot.');
    return;
end

% --- 2. Define window and extract data ---
if ~isfield(result, 'rate') || ~isfield(result.rate, 'mu') || isempty(result.rate.mu)
    warning('`result.rate.mu` is missing or empty. Cannot plot predicted firing rate.');
    return;
end
predicted_series = result.rate.mu;

start_sec = event_time_sec - options.Window(1);
duration_sec = sum(options.Window);

% --- 3. Prepare title and generate plot ---
plot_title = options.Title;
if strlength(plot_title) == 0 % Generate default title if none provided
    plot_title = sprintf('Model Prediction around #%d "%s" Event', ...
        options.EventIndex, options.EventType);
end

% Assuming `plot_predictions` is on the path and has the expected signature
h_fig = plot_predictions(result.stim, result.sps, predicted_series, result.events, ...
    start_sec, duration_sec, plot_title);

end