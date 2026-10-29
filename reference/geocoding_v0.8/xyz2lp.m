function [line,pixel] = xyz2lp(xyz_vec,image,orbit,...
                               MAXITER,CRITERPOS)
    
% $$$ % $$$ initialization
% $$$ num_points = size(xyz_vec,1);

% get times
[t_azi,t_ran] = xyz2t(xyz_vec,image,orbit,MAXITER,CRITERPOS);

% $$$ SLC.parameters
ta1 = image(1);    %second of day of first azimuth line
PRF = image(3);    %pulse repetition frequency
tr1 = image(2);    %one way time (s) to first range pixel
RSR = image(4);    %range sampling rate


line  = 1 + PRF.*(t_azi - ta1);
pixel = 1 + RSR.*(t_ran - tr1);

% $$$ EOF
