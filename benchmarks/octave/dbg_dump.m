args = argv(); toolbox = args{1}; scen = args{2}; outfile = args{3};
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'shim'));
addpath(fullfile(toolbox,'main')); addpath(fullfile(toolbox,'utils')); addpath(fullfile(toolbox,'lambda'));
ps_set_globals;
global project_id ps_method m2ph sat_vel althyp_index althyp_index2
project_id='dbg'; ps_method='boot'; m2ph=-4*pi/0.055517120823757024; sat_vel=7690;
althyp_index=NaN; althyp_index2=NaN;
S=load(fullfile(scen,'est_input.mat'));
[param,acheck,info,varfac,coh]=ps_ils_and_bootstrap_estimation(S.h2ph,S.Btemp,S.Bdop,S.dphase,S.model,S.std_param,S.sig2_est,S.breakpoint,S.breakpoint2);
dlmwrite(outfile,[varfac coh acheck]);

