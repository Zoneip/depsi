function [norm_orbit,norm_orbit_line] = intrp_orbit(line,image,precorb)

% $$$ intrp_orbit.m
% $$$  
% $$$ Made by Petar Marinkovic
% $$$ Login   <pmar@cake.lr.tudelft.nl>
% $$$  
% $$$ Started on  Thu Dec 01 12:00:00 2005 Petar Marinkovic
% $$$ Last update Wed Oct  3 14:40:55 2007 Petar Marinkovic
% $$$ 
% $$$ DESCRIPTION: normalizes and interpolates the input orbit. Function
% $$$ also dumps normalized XYZ point on the input orbit corresponding an
% $$$ input line (azimut ~ slow time). image variable represents metadata
% $$$ vector created with metadata.m.
% $$$ 
% $$$ NOTE.1: In interp1 functions 'LINEAR' interpolator used : with this one
% $$$ I obtained the most optimal results, in terms of residuals, if
% $$$ geocoding loop is done. Other options are splines of a certain
% $$$ degree. Strange?!?! 
% $$$ 
% $$$ NOTE.2: Function is suboptimal in cases where the large lines is used
% $$$ as input file, since for each line orbit is again interpolated and
% $$$ normalized.
% $$$ 
% $$$ INPUT: [line]   : line for which intrp.XYZ.orb is dumped
% $$$        [image]  : metadata vector
% $$$        [precorb]: orbit vector
% $$$ 
% $$$ OUTPUT: [norm_orbit]     : normalized interpolated orbit
% $$$         [norm_orbit_line]: normalized interpolated state vector 
% $$$                            for input line
% $$$
% $$$ TODO: ....
% $$$ 
% $$$ CHANGE.LOG:
% $$$ 
% $$$ [01.12.05, 12:00]: initial version and commit
% $$$ [13.03.07, 15:06]: support for input/output vectors

% interpolator
% $$$ intrpOrder = 'linear';
% $$$ intrpOrder = 'cubic';
intrpOrder = 'spline';

% compute normalization factors
px = precorb(:,1); % time
f  = min(px);
g  = (max(px)-min(px));
px = (px-f)/g;

% polyDegree
polyDegree = 2;

coef_x1 = (polyfit(px,precorb(:,2),polyDegree));
a = coef_x1(3);
b = coef_x1(2);
c = coef_x1(1);
coef_x = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_y1 = (polyfit(px,precorb(:,3),polyDegree));
a = coef_y1(3);
b = coef_y1(2);
c = coef_y1(1);
coef_y = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_z1 = (polyfit(px,precorb(:,4),polyDegree));
a = coef_z1(3);
b = coef_z1(2);
c = coef_z1(1);
coef_z = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

vel_x = polyval(polyder(coef_x),precorb(:,1));
vel_y = polyval(polyder(coef_y),precorb(:,1));
vel_z = polyval(polyder(coef_z),precorb(:,1));

acc_x = repmat(polyder(polyder(coef_x)),numel(precorb(:,1)),1);
acc_y = repmat(polyder(polyder(coef_y)),numel(precorb(:,1)),1);
acc_z = repmat(polyder(polyder(coef_z)),numel(precorb(:,1)),1);

% interpolated orbit
norm_orbit = [precorb(:,1),...
              precorb(:,2),precorb(:,3),precorb(:,4),...
              vel_x,  vel_y,  vel_z,...
              acc_x,  acc_y,  acc_z];

% interpolated orbit for l_aztime
l_aztime = (line-1)/image(3) + image(1);

% interpolate satellite position for l_aztime
pos = interp1(precorb(:,1),precorb(:,2:4),l_aztime,intrpOrder);
vel = interp1(precorb(:,1),[vel_x vel_y vel_z],l_aztime,intrpOrder);
acc = interp1(precorb(:,1),[acc_x acc_y acc_z],l_aztime,intrpOrder);

if length(line) == 1;
    norm_orbit_line = [l_aztime reshape(pos,1,numel(pos)) ...
                       reshape(vel,1,numel(vel)) ...
                       reshape(acc,1,numel(acc))];
else
    norm_orbit_line = [l_aztime pos vel acc];
end
    
% $$$ EOF
