args = argv(); toolbox = args{1}; scen = args{2};
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'shim'));
addpath(fullfile(toolbox,'main')); addpath(fullfile(toolbox,'utils')); addpath(fullfile(toolbox,'lambda'));
ps_set_globals;
global project_id ps_method m2ph sat_vel althyp_index althyp_index2
project_id='dbg'; ps_method='boot'; m2ph=-4*pi/0.055517120823757024; sat_vel=7690;
althyp_index=NaN; althyp_index2=NaN;
S=load(fullfile(scen,'est_input.mat')); Nifgs=size(S.dphase,2); Npar=2; model=1;
sig2=[(pi*15/180)^2; repmat((pi*20/180)^2,Nifgs,1)];
Qy1=repmat(2*sig2(1),Nifgs,Nifgs); for v=1:Nifgs, Qy1(v,v)=Qy1(v,v)+2*sig2(v+1); end
printf('DBG_SIG2 %.10f %.10f\n', sig2(1), sig2(2));
printf('DBG_QY1DIAG %.10f %.10f\n', Qy1(1,1), Qy1(2,2));
invQy1=inv(Qy1);
mean_h2ph=mean(S.h2ph,1)';
[B1Qy2,par_index]=ps_model_definitions('model',model,Nifgs,mean_h2ph,S.Btemp,S.Bdop,S.std_param);
B1=B1Qy2(1:Nifgs,:);
rhs_fixed=inv(B1'*invQy1*B1)*B1'*invQy1;
y=S.dphase(1,:)'; bcheck=rhs_fixed*y; e=y-B1*bcheck;
printf('DBG_E '); printf('%.8f ', e(1:5)); printf('\n');
printf('DBG_OMT %.10f\n', e'*invQy1*e);
printf('DBG_RHS_ROW1 '); printf('%.6e ', rhs_fixed(1,1:4)); printf('\n');
printf('DBG_B1COL1 '); printf('%.6f ', B1(1:5,1)'); printf('\n');
