testPath = '/Users/matt/Documents/Zhao Lab/audio/test clips';

wavPath = fullfile(testPath, 'S178_test.wav');
heardPath = fullfile(testPath, 'S178_test_heard.txt');
producedPath = fullfile(testPath, 'S178_test_produced.txt');

outputPath = fullfile(testPath, 'S178_test_state.txt');

export_state_labels_to_audacity_function(wavPath, heardPath, producedPath, outputPath);

function export_state_labels_to_audacity_function(wavPath, heardLabelPath, producedLabelPath, outputPath, configPath)
% export_state_labels_to_audacity emit conversational states as audacity labels.
%
% export_state_labels_to_audacity(wavPath, heardLabelPath, producedLabelPath, outputPath)
% loads a wav file alongside separate audacity txt label files for perceived
% (heard) and produced events. it computes conversational/spontaneous states
% using project helpers and writes an audacity-compatible label track with
% "convo" and "spon" segments. an optional config json path can be supplied
% as a fifth argument; defaults.json is used otherwise.

%% locate repository root and add required paths
% we ensure src utilities are available before invoking helpers.
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
addpath(genpath(fullfile(repoRoot, 'src')));

%% validate inputs
% we require four positional arguments and normalize optional config.
if nargin < 4
    error('glm:InvalidInput', ['Usage: export_state_labels_to_audacity(wavPath, heardLabelPath, ', ...
        'producedLabelPath, outputPath, [configPath])']);
end
if nargin < 5 || isempty(configPath)
    configPath = fullfile(repoRoot, 'config', 'defaults.json');
end

wavPath = ensure_char_path(wavPath, 'wavPath');
heardLabelPath = ensure_char_path(heardLabelPath, 'heardLabelPath');
producedLabelPath = ensure_char_path(producedLabelPath, 'producedLabelPath');
outputPath = ensure_char_path(outputPath, 'outputPath');
configPath = ensure_char_path(configPath, 'configPath');

if exist(wavPath, 'file') ~= 2
    error('glm:FileNotFound', 'Audio file not found: %s', wavPath);
end
if exist(heardLabelPath, 'file') ~= 2
    error('glm:FileNotFound', 'Heard label file not found: %s', heardLabelPath);
end
if exist(producedLabelPath, 'file') ~= 2
    error('glm:FileNotFound', 'Produced label file not found: %s', producedLabelPath);
end
if exist(configPath, 'file') ~= 2
    error('glm:FileNotFound', 'Config file not found: %s', configPath);
end

%% load configuration and supporting data
% we decode defaults, grab state parameters, load label structs, and fetch audio duration.
config = jsondecode(fileread(configPath));
if ~isfield(config, 'dt')
    error('glm:InvalidConfig', 'Config must define a dt field.');
end
if ~isfield(config, 'state')
    error('glm:InvalidConfig', 'Config must define a state struct.');
end
stateCfg = config.state;
dt = double(config.dt);
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('glm:InvalidConfig', 'Config dt must be a positive finite scalar.');
end

audioInfo = audioinfo(wavPath);
audioDuration = audioInfo.Duration;

heardEvents = load_labels(heardLabelPath, 'perceived');
producedEvents = load_labels(producedLabelPath, 'produced');
allEvents = [heardEvents(:); producedEvents(:)];

%% build stimulus grid covering audio duration and events
% we create a uniform grid spanning the greater of audio length and latest label offset.
maxEventTime = infer_max_event_time(allEvents);
T = max(audioDuration, maxEventTime);
numBins = floor(T / dt + 0.5) + 1;
t = (0:numBins-1)' * dt;
stim = struct('t', t, 'dt', dt, 'mask', struct('good', true(numBins, 1)));

%% compute conversational states on the grid
% we reuse project helper logic to derive convo and spon boolean streams.
states = compute_states(allEvents, stim, stateCfg);

%% convert states to audacity intervals
% we stitch contiguous bins of the same state into labelled intervals ready for export.
segments = build_state_segments(states, stim);

%% write the audacity label file
% we export tab-delimited rows representing each labelled state segment.
write_audacity_labels(outputPath, segments);

fprintf('wrote %d state segments to %s\n', numel(segments), outputPath);
end

function pathOut = ensure_char_path(pathIn, fieldName)
% coerce supported string types into a char vector and validate input.
if isstring(pathIn)
    if numel(pathIn) ~= 1
        error('glm:InvalidInput', '%s must be a single string scalar.', fieldName);
    end
    pathOut = char(pathIn);
elseif ischar(pathIn)
    pathOut = pathIn;
else
    error('glm:InvalidInput', '%s must be a character vector or string scalar.', fieldName);
end
end

function maxEventTime = infer_max_event_time(events)
% derive the latest event offset to ensure the grid covers all annotations.
if isempty(events)
    maxEventTime = 0;
    return
end
if ~isstruct(events) || ~isfield(events, 't_off')
    error('glm:InvalidInput', 'Events must contain t_off fields.');
end
values = [events.t_off];
if isempty(values)
    maxEventTime = 0;
else
    maxEventTime = max(double(values));
end
end

function segments = build_state_segments(states, stim)
% convert boolean convo/spon vectors into contiguous labelled intervals.

%% prepare per-bin labels restricted to valid bins
% we assign textual labels for conversational and spontaneous bins.
convo = logical(states.convo(:));
spon = logical(states.spon(:));
good = convo | spon;
labels = strings(numel(convo), 1);
labels(convo) = "convo";
labels(spon) = "spon";

%% iterate through bins to coalesce contiguous runs of identical state
% we gather start/end indices for each labelled segment.
segments = repmat(struct('start', 0, 'stop', 0, 'label', ""), 0, 1);
if ~any(good)
    return
end

currentLabel = "";
startIdx = -1;
for ii = 1:numel(labels)
    label = labels(ii);
    if label == ""
        if startIdx ~= -1
            segments(end+1) = finalize_segment(stim, startIdx, ii-1, currentLabel); %#ok<AGROW>
            startIdx = -1;
            currentLabel = "";
        end
        continue
    end

    if startIdx == -1
        startIdx = ii;
        currentLabel = label;
    elseif label ~= currentLabel
        segments(end+1) = finalize_segment(stim, startIdx, ii-1, currentLabel); %#ok<AGROW>
        startIdx = ii;
        currentLabel = label;
    end
end

if startIdx ~= -1
    segments(end+1) = finalize_segment(stim, startIdx, numel(labels), currentLabel); %#ok<AGROW>
end
end

function segment = finalize_segment(stim, startIdx, endIdx, label)
% convert bin indices to time boundaries with the associated label.
dt = stim.dt;
segment = struct( ...
    'start', stim.t(startIdx), ...
    'stop', stim.t(endIdx) + dt, ...
    'label', label ...
);
end

function write_audacity_labels(outputPath, segments)
% emit audacity label file with conversational/spontaneous intervals.
fid = fopen(outputPath, 'w');
if fid < 0
    error('glm:FileIO', 'Unable to open output file: %s', outputPath);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

if isempty(segments)
    warning('glm:NoSegments', 'No labelled segments produced; empty file written.');
    return
end

for ii = 1:numel(segments)
    fprintf(fid, '%.6f\t%.6f\t%s\n', segments(ii).start, segments(ii).stop, char(segments(ii).label));
end
end
