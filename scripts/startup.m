% this script prepares the matlab path and required output folders.

repo_root = fileparts(fileparts(mfilename('fullpath')));
src_dir = fullfile(repo_root, 'src');

% load the src tree onto the matlab path for package resolution.
if ~isfolder(src_dir)
    error('vglm:startup:MissingSrcDir', 'expected src directory at %s', src_dir);
end
addpath(genpath(src_dir));

% make sure results and figures folders exist for downstream artifacts.
required_dirs = {'results', 'figures'};
for k = 1:numel(required_dirs)
    target_dir = fullfile(repo_root, required_dirs{k});
    if ~isfolder(target_dir)
        mkdir(target_dir);
    end
    assert(isfolder(target_dir), 'vglm:startup:MissingDir', ...
        'expected directory %s to be present after startup', target_dir);
end
