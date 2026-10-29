% bench_selection.m -- time ps_selection on a full DePSI-format project
%
% Usage: octave --no-gui --quiet bench_selection.m <toolbox_dir> <scenario_dir> <workdir>

args = argv();
if numel(args) < 3
  error('usage: bench_selection.m <toolbox_dir> <scenario_dir> <workdir>');
end
toolbox = args{1};
scen = args{2};
workdir = args{3};
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'shim'));

pkg load image;
pkg load statistics;

addpath(fullfile(toolbox, 'main'));
addpath(fullfile(toolbox, 'utils'));
addpath(fullfile(toolbox, 'lambda'));
addpath(fullfile(toolbox, 'quadtree'));

S = load(fullfile(scen, 'sel_input.mat'));

% ps_set_globals declares the globals; any of them that is *also* a function
% argument would otherwise shadow the argument with an empty global, so they all
% have to be filled in here (depsi.m does that via ps_readinput_parameters).
ps_set_globals;
global fig ps_method ps_eval_method project_id results_id
global m2ph Npar_max unwrap_istep unwrap2_istep
global max_mem_buffer orbit Nlines Npixels az_spacing r_spacing
global Npsc_selections weighting varfac_threshold
global visible_plots detail_plots ens_coh_threshold
global atmo_range defo_range sat_vel weighted_unwrap
global detrend_flag defo_model_flag run_mode
Npsc_selections = S.Npsc_selections;
weighting = 'vce'; varfac_threshold = 3; ens_coh_threshold = 0.9;
weighted_unwrap = 'yes'; detrend_flag = 'no'; defo_model_flag = 'no';
run_mode = 'normal'; results_id = 'orig'; orbit = 'asc';
atmo_range = 5000; defo_range = 5000; Npar_max = 10;
unwrap_istep = 1; unwrap2_istep = 1;
project_id = 'bench';
az_spacing = S.az_spacing;
r_spacing = S.r_spacing;
Nlines = S.Nlines;
Npixels = S.Npixels;
max_mem_buffer = S.max_mem_buffer;
ps_eval_method = 'psp';
visible_plots = 'n';
detail_plots = 'n';
fig = 0;
global m2ph sat_vel
m2ph = -4*pi/0.055517120823757024;
sat_vel = 7690;

if ~exist(workdir, 'dir'); mkdir(workdir); end
old = pwd; cd(workdir);

tic;
[grid_array_az, grid_array_r, Npsc, Npsp] = ps_selection( ...
    S.slc, S.ifg, S.h2ph, ...
    S.psc_selection_method, S.psc_selection_gridsize, S.psc_threshold, ...
    S.Npsc_selections, S.psp_selection_method, S.psp_threshold1, [], ...
    S.slc_selection, S.Btemp, [], [], S.crop, S.crop_in, ...
    'yes', S.calfactors, S.gamma_threshold, S.psc_distribution, ...
    S.livetime_threshold, S.peak_tolerance, [NaN NaN]);
t = toc;

cd(old);
printf('SEL_TIME %.4f\n', t);
printf('SEL_NPSC %d\n', Npsc(1));
printf('SEL_NPSP %d\n', Npsp(1));

% peak resident set size of this Octave process (no /usr/bin/time on this host)
s = fileread('/proc/self/status');
tok = regexp(s, 'VmHWM:\s+(\d+) kB', 'tokens');
if ~isempty(tok)
  printf('PEAK_RSS_MB %.1f\n', str2double(tok{1}{1})/1024);
end
