function [t_azi,t_ran] = xyz2t(xyz_vec,image,orbit,MAXITER,CRITERPOS)

np    = size(xyz_vec,1);

t_azi = zeros(np,1);
t_ran = zeros(np,1);

SOL = 299792458;   % speed of light

ta1         = image(1); % second of day of first azimuth line
PRF         = image(3); % pulse repetition frequency
midCropLine = image(16);% mid line as start point for interpolation

% $$$ intrpOrder  = 'linear';
% $$$ intrpOrder  = 'cubic';
intrpOrder  = 'spline';


% interpolate and normalize orbit for [midLine]
[orbit,orbit_line] = intrp_orbit(midCropLine,image,orbit);

tazi = (midCropLine-1)/PRF+ta1;

for k = 1 : np;
    
% $$$     if mod(k,50)==0
% $$$         disp(['Point ' num2str(k) ' of ' num2str(np) ' points in total'])
% $$$     end
    
    pos = xyz_vec(k,:);
    
    for iter = 1 : MAXITER;
        
        % update equations
        possat(1) = interp1(orbit(:,1),orbit(:,2),tazi);  %orbit corresponding to the azimuth time
        possat(2) = interp1(orbit(:,1),orbit(:,3),tazi);
        possat(3) = interp1(orbit(:,1),orbit(:,4),tazi);
        velsat(1) = interp1(orbit(:,1),orbit(:,5),tazi);  
        velsat(2) = interp1(orbit(:,1),orbit(:,6),tazi);
        velsat(3) = interp1(orbit(:,1),orbit(:,7),tazi);
        accsat(1) = interp1(orbit(:,1),orbit(:,8),tazi);  
        accsat(2) = interp1(orbit(:,1),orbit(:,9),tazi);
        accsat(3) = interp1(orbit(:,1),orbit(:,10),tazi);

        delta = pos - possat;
        
        % update solution
        solution = -(velsat(1)*delta(1)+velsat(2)*delta(2)+velsat(3)* ...
                     delta(3))/(accsat(1)*delta(1)+accsat(2)*delta(2)+ ...
                                accsat(3)*delta(3)-velsat(1)^2-velsat(2)^2- ...
                                velsat(3)^2);

        tazi = tazi + solution;
        
        % check convergence
        if(abs(solution) < CRITERPOS)
            break;
        else
            if ( iter >= MAXITER )
                MAXITER=MAXITER+1;
            end
        end
        
    end

    t_azi(k) = tazi; % azimuth time
    
    %orbit corresponding to the azimuth time
    possat(1) = interp1(orbit(:,1),orbit(:,2),t_azi(k),intrpOrder); 
    possat(2) = interp1(orbit(:,1),orbit(:,3),t_azi(k),intrpOrder);
    possat(3) = interp1(orbit(:,1),orbit(:,4),t_azi(k),intrpOrder);
    
    delta = pos - possat;
    
    t_ran(k) = sqrt(delta(1)^2+delta(2)^2+delta(3)^2)/SOL;  %range time, tran
    
end

% EOF
