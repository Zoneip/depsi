function [xyz] = lph2xyz(line,pixel,height,ellipsoid,image,norm_orbit_line, ...
                         MAXITER,CRITERPOS)

% $$$ initialization
num_points = numel(line);
xyz = zeros(num_points,3);

% $$$ parmeters of the image
tr1          = image(2); % one.way range time [sec]
RSR          = image(4); % one.way in [Hz]
centerphi    = image(6); 
centerlambda = image(7);

SOL = 299792458;

% $$$ reference surface: WGS84

ell_a  = ellipsoid(1);                 % semimajor of the ellipsoid
ell_b  = ellipsoid(2);                 % semiminor of the ellipsoid
ell_e2 = (ell_a^2-ell_b^2)/ell_a^2;   % squared first eccentricity(derived)
% $$$ ell_e2b=(ell_a^2-ell_b^2)/ell_b^2;  % squared second eccentricity(derived)

% $$$ [lat,long,h] of scene center to [x,y,z]
h            = mean(height); % this is only for initial values
centerphi    = centerphi*pi/180;
centerlambda = centerlambda*pi/180;
Ncenter      = ell_a/sqrt(1-ell_e2*(sin(centerphi)^2)); % eccentricity
scenecenterx = (Ncenter+h)*cos(centerphi)*cos(centerlambda);
scenecentery = (Ncenter+h)*cos(centerphi)*sin(centerlambda);
scenecenterz = (Ncenter+h-ell_e2*Ncenter)*sin(centerphi);

for n = 1 : num_points; % loop through points

    posonellx = scenecenterx;
    posonelly = scenecentery;
    posonellz = scenecenterz;
    
    % convert line and pixel to time: from Doris
    % aztime = line2ta(line,image.t_zai1,image.prf)
    % ratime = pix2tr(pixel,image.t_range1,image.rsr2x)
    % $$$     aztime = ta1 + (line(n) -1.0)/PRF; % not used?
    ratime = tr1 + (pixel(n)-1.0)/RSR;
    
    %get position and velocity of the satellite
    possatx = norm_orbit_line(n,2);
    possaty = norm_orbit_line(n,3);
    possatz = norm_orbit_line(n,4);
    velsatx = norm_orbit_line(n,5);
    velsaty = norm_orbit_line(n,6);
    velsatz = norm_orbit_line(n,7);
    
    for iter=1:MAXITER
        
        %update equations and slove system
        dsat_Px = posonellx - possatx;   %vector of 'satellite to P on ellipsoid'
        dsat_Py = posonelly - possaty;
        dsat_Pz = posonellz - possatz;
        
        equationset(1,1) = -(velsatx*dsat_Px+velsaty*dsat_Py+velsatz* ...
                             dsat_Pz);
        
        equationset(2,1) = -(dsat_Px*dsat_Px+dsat_Py*dsat_Py+dsat_Pz* ...
                             dsat_Pz-(SOL*ratime).^2);
        
        equationset(3,1) = -((posonellx*posonellx+posonelly*posonelly)/ ...
                             ((ell_a+height(n)).^2)+(posonellz/(ell_b+ ...
                                                          height(n)))^2-1.0);

        partialsxyz(1,1) = velsatx;
        partialsxyz(1,2) = velsaty;
        partialsxyz(1,3) = velsatz;
        partialsxyz(2,1) = 2*dsat_Px;
        partialsxyz(2,2) = 2*dsat_Py;
        partialsxyz(2,3) = 2*dsat_Pz;
        partialsxyz(3,1) = (2*posonellx)/((ell_a+height(n))^2);
        partialsxyz(3,2) = (2*posonelly)/((ell_a+height(n))^2);
        partialsxyz(3,3) = (2*posonellz)/((ell_b+height(n))^2);
        
        % solve system [NOTE] orbit_normalized, otherwise close to
        % singular
        solpos = partialsxyz\equationset;
        solx = solpos(1);
        soly = solpos(2);
        solz = solpos(3);
        
        % update solution
        posonellx = posonellx + solx;
        posonelly = posonelly + soly;
        posonellz = posonellz + solz;
        
        % check convergence
        if (abs(solx)<CRITERPOS && abs(soly)<CRITERPOS && abs(solz)<CRITERPOS)
            break;
        else
            if(iter>=MAXITER)
                MAXITER=MAXITER+1;
            end
        end
        
    end
    
    % final solution: array of XYZ coordinates
    xyz(n,:)=[posonellx posonelly posonellz];
    
end

% $$$ EOF
