function [t_azi,t_ran] = lp2t(line,pixel,image)

% assumes [line,pixel] are in "absolute image" coordinates, not cropped
    
% $$$ parmeters of the image
ta1          = image(1); % start time in UTC [sec]
tr1          = image(2); % one.way range time [sec]
PRF          = image(3); % [Hz]
RSR          = image(4); % one.way in [Hz]

t_azi = ta1 + (line - 1.0)./PRF; % not used?
t_ran = tr1 + (pixel - 1.0)./RSR;
    
% $$$ EOF
