function capture_baseline(config_file)
% Native MATLAB only. Run in a fresh process with matlab -batch.
% No Octave compatibility functions are added. Production source is unchanged.
% A completed entry is NOT proof that every requested processing branch ran.
cfg = jsondecode(fileread(config_file));
if exist('OCTAVE_VERSION', 'builtin')
    error('depsipy:NativeMATLABRequired', 'Use native MATLAB, not the Octave shim route.');
end
old_path = path;
old_dir = pwd;
cleanup = onCleanup(@() restore_session(old_path, old_dir)); %#ok<NASGU>
% Remove any inherited path into our Octave shim directory.
paths = strsplit(path, pathsep);
for k = 1:numel(paths)
    normalized = strrep(paths{k}, '\', '/');
    if contains(normalized, '/benchmarks/octave/shim')
        rmpath(paths{k});
    end
end
[entry_dir, entry_name] = fileparts(cfg.entry_file);
root = fileparts(entry_dir); % supplied entry convention: <toolbox>/main/depsi.m
for name = {'utils','lambda','quadtree','simulation','geocoding_v0.8'}
    p = fullfile(root, name{1});
    if isfolder(p), addpath(p); end
end
addpath(entry_dir, '-begin');
resolved = which(entry_name);
if ~strcmp(resolved, cfg.entry_file)
    error('depsipy:WrongEntry', 'Resolved entry differs from config: %s', resolved);
end
cd(cfg.work_dir);
report = struct('schema', 1, 'status', 'started', 'matlab_version', version, ...
    'release', version('-release'), 'computer', computer, 'toolboxes', ver, ...
    'entry', resolved, 'function_resolution', struct(), 'profile_enabled', cfg.profile);
for name = {'gcp','parpool','feature','nanmean','nanstd','ps_densification','ps_noise_filter','ps_spatio_temporal_consistency'}
    report.function_resolution.(name{1}) = which(name{1});
end
report.parallel_before = pool_information();
report.started_utc = char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ssXXX"));
write_json(cfg.result_file, report);
if cfg.profile, profile clear; profile on; end
start_time = tic;
try
    % Function entry accepts an explicit parameter-file path. Script wrappers
    % with hardcoded work_dir require a separately reviewed production adapter.
    feval(entry_name, cfg.param_file);
    report.status = 'completed';
catch err
    report.status = 'failed';
    report.error = struct('identifier', err.identifier, 'message', err.message, ...
        'report', getReport(err, 'extended', 'hyperlinks', 'off'));
end
report.elapsed_s = toc(start_time);
report.parallel_after = pool_information();
if cfg.profile
    profile off;
    profile_info = profile('info');
    save(cfg.profile_file, 'profile_info');
    names = {profile_info.FunctionTable.FunctionName};
    report.profile_function_names = names;
    % Call counts are evidence, not claims of scientific success.
    report.calls = struct();
    for name = {'ps_network','ps_ils_and_bootstrap_estimation','ps_densification','ps_noise_filter','ps_spatio_temporal_consistency'}
        match = contains(names, name{1});
        report.calls.(name{1}) = sum([profile_info.FunctionTable(match).NumCalls]);
    end
else
    report.calls_unverified = 'Enable --profile for function call evidence; outer entry completion is insufficient.';
end
report.artifacts = dir(fullfile(cfg.work_dir, '**', '*'));
write_json(cfg.result_file, report);
if strcmp(report.status, 'failed')
    error('depsipy:BaselineFailed', '%s', report.error.message);
end
end

function data = pool_information()
data = struct('available', false, 'workers', 0);
if exist('gcp', 'file')
    data.available = true;
    try
        pool = gcp('nocreate');
        if ~isempty(pool), data.workers = pool.NumWorkers; end
    catch err
        data.error = err.message;
    end
end
end

function write_json(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('depsipy:OutputFailure', 'Cannot write %s', filename); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(value));
end

function restore_session(old_path, old_dir)
path(old_path);
cd(old_dir);
end
