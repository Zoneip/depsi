function [x,y,z]=ell2xyz(ell,phi,lambda,height)
% convert geodetic coordinate phi,lambda,height to geocentered cartesian
% coordinate
%
% output:
% x,y,z              ---  geocentered cartesian coordinate
%
% input:
% ell(1),ell(2)      ---  semimajor amd semiminor axis of ellipsoid
% phi,lambda,height  ---  geodetic coordinate,latitude,longitude(in degree),height
%********************************************

if (nargin<1); error('ell2xyz: nargin=0'); end;

% $$$ WGS84: ellipsoid
% $$$ a=	6378137.000
% $$$ b=	6356752.314
% $$$ e2=	0.00669437999
% $$$ 1/f=	298.2572235

ell_a = ell(1);
ell_b = ell(2);

phi    =  (phi/180)*pi;    % [radians]
lambda = (lambda/180)*pi;  % [radians]

ell_e2 = (ell_a^2-ell_b^2)/ell_a^2; % squared first eccentricity

N   = ell_a ./ sqrt(1-ell_e2.*(sin(phi)).^2);
Nph = N + height;
x   =  Nph .* cos(phi)  .* cos(lambda); 
y   =  Nph .* cos(phi)  .* sin(lambda); 
z   = (Nph - ell_e2*N) .* sin(phi);

% EOF
