args = argv(); toolbox = args{1}; scen = args{2}; arc = str2double(args{3});
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'shim'));
addpath(fullfile(toolbox,'main')); addpath(fullfile(toolbox,'utils')); addpath(fullfile(toolbox,'lambda'));
ps_set_globals;
global project_id ps_method m2ph sat_vel althyp_index althyp_index2
project_id='dbg'; ps_method='boot'; m2ph=-4*pi/0.055517120823757024; sat_vel=7690;
althyp_index=NaN; althyp_index2=NaN;
S=load(fullfile(scen,'est_input.mat')); Nifgs=size(S.dphase,2); Npar=2; model=1;
sig2=[(pi*15/180)^2; repmat((pi*20/180)^2,Nifgs,1)];
Qy1=repmat(2*sig2(1),Nifgs,Nifgs); for v=1:Nifgs, Qy1(v,v)=Qy1(v,v)+2*sig2(v+1); end
mean_h2ph=mean(S.h2ph,1)';
[B1Qy2,par_index]=ps_model_definitions('model',model,Nifgs,mean_h2ph,S.Btemp,S.Bdop,S.std_param);
B1=B1Qy2(1:Nifgs,:); Qy2=B1Qy2(Nifgs+1:end,:);
A=[-2*pi*eye(Nifgs); zeros(Npar,Nifgs)]; B=[B1;eye(Npar)]; C=[A B];
Qy=[Qy1 zeros(Nifgs,Npar); zeros(Npar,Nifgs) Qy2];
Qahat=inv(C'*inv(Qy)*C); Qahat=Qahat(1:Nifgs,1:Nifgs);
[Z,L,D]=decorrel_freek(Qahat);
afloat=S.dphase(arc+1,:)'/(-2*pi); zfloat=Z'*afloat;
[Chi2,zboot]=ps_chistart(D,L,zfloat,1);
printf('ARC_CHI2 %.12f\n', Chi2);
printf('ARC_ZBOOT '); printf('%.0f ', zboot); printf('\n');
printf('ARC_ZFLOAT '); printf('%.9f ', zfloat(1:6)); printf('\n');
printf('ARC_D '); printf('%.9f ', D(1:6)); printf('\n');
