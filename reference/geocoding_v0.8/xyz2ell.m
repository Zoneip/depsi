function [phi_lam_height] = xyz2ell(position,ellipsoid)

ell_a   = ellipsoid(1);
ell_b   = ellipsoid(2);

ell_e2  = (ell_a^2-ell_b^2)/ell_a^2;    % squared first eccentricity(derived)
ell_e2b = (ell_a^2-ell_b^2)/ell_b^2;    % squared second eccentricity(derived)

posx = position(:,1);
posy = position(:,2);
posz = position(:,3);

r    = sqrt(posx.^2 + posy.^2);
mu   = atan2(posz.*ell_a, r.*ell_b);

sin3 = (sin(mu)).^3;
cos3 = (cos(mu)).^3;
phi  = atan2((posz + ell_e2b .* ell_b .* sin3),...
             (r - ell_e2 .* ell_a .* cos3));

lambda = atan2(posy,posx);
N = ell_a ./ sqrt(1 - ell_e2 .* (sin(phi)).^2); % for every point no
                                                % approx with scene.center
height = (r ./ cos(phi)) - N;

phi_lam_height = [(phi.*180)./pi (lambda.*180)./pi height];
    
% $$$ EOF
