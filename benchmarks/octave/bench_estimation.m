% bench_estimation.m -- time the per-arc ambiguity resolution (LAMBDA/ILS)
%
% Usage: octave --no-gui --quiet bench_estimation.m <toolbox_dir> <scenario_dir>
% Prints machine-readable lines "KEY value".

args = argv();
if numel(args) < 2
  error('usage: bench_estimation.m <toolbox_dir> <scenario_dir>');
end
toolbox = args{1};
scen = args{2};
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'shim'));

addpath(fullfile(toolbox, 'main'));
addpath(fullfile(toolbox, 'utils'));
addpath(fullfile(toolbox, 'lambda'));
addpath(fullfile(toolbox, 'quadtree'));

ps_set_globals;
global project_id ps_method ens_coh_threshold varfac_threshold althyp_index althyp_index2
project_id = 'bench';
ps_method = 'boot';
ens_coh_threshold = 0.9;
varfac_threshold = 3;
althyp_index = NaN;
althyp_index2 = NaN;
global m2ph sat_vel
m2ph = -4*pi/0.055517120823757024;
sat_vel = 7690;


S = load(fullfile(scen, 'est_input.mat'));
printf('EST_ARCS %d\n', size(S.dphase, 1));
printf('EST_NIFGS %d\n', size(S.dphase, 2));

% warm-up (JIT / function loading)
ps_ils_and_bootstrap_estimation(S.h2ph(1:20,:), S.Btemp, S.Bdop, S.dphase(1:20,:), ...
                                S.model, S.std_param, S.sig2_est, ...
                                S.breakpoint, S.breakpoint2);

tic;
[param, acheck, info, varfac, coh] = ps_ils_and_bootstrap_estimation( ...
    S.h2ph, S.Btemp, S.Bdop, S.dphase, S.model, S.std_param, S.sig2_est, ...
    S.breakpoint, S.breakpoint2);
t = toc;

printf('EST_TIME %.4f\n', t);
printf('EST_FINITE_VARFAC %d\n', sum(isfinite(varfac)));
printf('EST_MEDIAN_VARFAC %.4f\n', median(varfac(isfinite(varfac))));
printf('EST_MEDIAN_COH %.4f\n', median(coh(isfinite(coh))));

% peak resident set size of this Octave process (no /usr/bin/time on this host)
s = fileread('/proc/self/status');
tok = regexp(s, 'VmHWM:\s+(\d+) kB', 'tokens');
if ~isempty(tok)
  printf('PEAK_RSS_MB %.1f\n', str2double(tok{1}{1})/1024);
end
