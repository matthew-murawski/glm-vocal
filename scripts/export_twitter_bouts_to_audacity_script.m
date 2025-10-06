samplePath = '/Users/matt/Documents/GitHub/vocalization/data/Label Files/S177';

producedLabelPath = fullfile(samplePath, 'M93A_S177_produced.txt');
outputLabelPath = fullfile(samplePath, 'M93A_S177_twitter.txt');

export_twitter_bouts_to_audacity(producedLabelPath, outputLabelPath);

function export_twitter_bouts_to_audacity(producedPath, outputPath, configPath)
% export_twitter_bouts_to_audacity collapse twitter syllables and write labels.
%
% export_twitter_bouts_to_audacity(producedPath, outputPath) loads a
% produced-label audacity txt file, collapses produced twitter syllables
% into bouts using consolidate_twitter_bouts, and emits a new audacity label
% file containing only the consolidated twitter bouts. an optional config
% json path can be supplied to override the default bout window.

%% locate repository root and add source helpers
% we ensure project src utilities are discoverable before invoking helpers.
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
addpath(genpath(fullfile(repoRoot, 'src')));

%% validate required arguments and normalize optional inputs
% we require produced and output paths; config defaults to config/defaults.json.
if nargin < 2
    error('glm:InvalidInput', 'Usage: export_twitter_bouts_to_audacity(producedPath, outputPath, [configPath])');
end
if nargin < 3 || isempty(configPath)
    configPath = fullfile(repoRoot, 'config', 'defaults.json');
end

producedPath = ensure_char_path(producedPath, 'producedPath');
outputPath = ensure_char_path(outputPath, 'outputPath');
configPath = ensure_char_path(configPath, 'configPath');

if exist(producedPath, 'file') ~= 2
    error('glm:FileNotFound', 'Produced label file not found: %s', producedPath);
end
if exist(configPath, 'file') ~= 2
    error('glm:FileNotFound', 'Config file not found: %s', configPath);
end

%% load configuration and extract twitter bout window
% we decode the defaults json and ensure a valid positive window is provided.
config = jsondecode(fileread(configPath));
if ~isfield(config, 'twitter_bout_window_s')
    error('glm:InvalidConfig', 'Config must define twitter_bout_window_s.');
end
boutWindow = double(config.twitter_bout_window_s);
if ~isscalar(boutWindow) || ~isfinite(boutWindow) || boutWindow <= 0
    error('glm:InvalidConfig', 'twitter_bout_window_s must be a positive finite scalar.');
end

%% load produced events and consolidate twitter syllables
% we parse the audacity txt file, collapse syllables into bouts, and retain produced twitter events only.
producedEvents = load_labels(producedPath, 'produced');
consolidated = consolidate_twitter_bouts(producedEvents, boutWindow);

if isempty(consolidated)
    boutEvents = consolidated;
else
    isProduced = strcmpi({consolidated.kind}, 'produced');
    labelStrings = string({consolidated.label});
    labelStrings = lower(strtrim(labelStrings(:)));
    labelStrings = replace(labelStrings, "_", "");
    labelStrings = replace(labelStrings, "-", "");
    labelStrings = replace(labelStrings, " ", "");
    isTwitter = contains(labelStrings, "twitter");

    boutEvents = consolidated(isProduced(:) & isTwitter(:));
end

%% convert bout events into audacity segments
% we transform each consolidated bout into a labelled interval ready for export.
segments = repmat(struct('start', 0, 'stop', 0, 'label', ""), 0, 1);
for ii = 1:numel(boutEvents)
    event = boutEvents(ii);
    tOn = double(event.t_on);
    tOff = double(event.t_off);
    if isempty(tOff) || ~isfinite(tOff) || tOff < tOn
        tOff = tOn;
    end

    segments(end+1, 1) = struct( ...
        'start', tOn, ...
        'stop', tOff, ...
        'label', "twitter_bout" ...
    ); %#ok<AGROW>
end

%% write the audacity label file
% we persist the consolidated bouts to a tab-delimited label track.
write_audacity_segments(outputPath, segments);

fprintf('wrote %d twitter bouts to %s\n', numel(segments), outputPath);
end

function pathOut = ensure_char_path(pathIn, fieldName)
% coerce supported string types into char vectors and guard against invalid inputs.
%% normalize supported input types
% we coerce strings and chars to a consistent representation and enforce scalars.
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

function write_audacity_segments(outputPath, segments)
% emit a tab-delimited audacity label file for the provided segments.
%% open the output file safely
% we ensure the destination is writable before streaming segments.
fid = fopen(outputPath, 'w');
if fid < 0
    error('glm:FileIO', 'Unable to open output file: %s', outputPath);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

%% write each segment as an audacity label row
% we warn when no segments exist then serialize one row per bout.
if isempty(segments)
    warning('glm:NoSegments', 'No twitter bouts found; writing empty label file.');
    return
end

for ii = 1:numel(segments)
    fprintf(fid, '%.6f\t%.6f\t%s\n', segments(ii).start, segments(ii).stop, char(segments(ii).label));
end
end
