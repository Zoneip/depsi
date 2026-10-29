function [container,orbit,containerChar] = metadataNew(result_file)

% $$$ metadata.m
% $$$  
% $$$ Made by Petar Marinkovic
% $$$ Login   <pmar@cake.lr.tudelft.nl>
% $$$  
% $$$ Started on  Thu Dec 01 12:00:00 2005 Petar Marinkovic
% $$$ Last update Wed Oct  3 14:36:12 2007 Petar Marinkovic
% $$$ 
% $$$ DESCRIPTION: normalizes and interpolates the input orbit. Function
% $$$ also dumps normalized XYZ point on the input orbit corresponding an
% $$$ input line (azimut ~ slow time). container variable represents metadata
% $$$ vector created with metadata.m.
% $$$ 
% $$$ INPUT: [result_file]: filename (including path) of Doris result file
% $$$ 
% $$$ TODO: ....
% $$$ 
% $$$ CHANGE.LOG:
% $$$ 
% $$$ [01.12.05, 12:00]: initial version and commit
% $$$ [12.03.07, 19:15]: flags and cheks for oversampling
% $$$ [01.04.08, 19:20]: new variable for sensor wavelength
% $$$ [21.04.08, 14:38]: support for ALOS res files
% $$$ [21.04.08, 16:38]: doppler info listed
% $$$ [23.04.08, 14:36]: different CONTAINER vector structure: cells

% $$$ CONTAINER.vector structure ---------------------------------
% $$$ --------------------------------------------------------
% $$$ container(1): acquisition start time
% $$$ container(2): 1way range time to 1st pix in [sec]
% $$$ container(3): PRF [Hz]
% $$$ container(4): RSR [Hz]: 2way!
% $$$ container(5): sensor wavelength [m]
% $$$ container(6): center pixel latitude
% $$$ container(7): center pixel longitude
% $$$ container(8): ovs.factor: AZIMUTH
% $$$ container(9): ovs.factor: RANGE
% $$$ container(10): Number_of_lines_original
% $$$ container(11): Number_of_pixels_original
% $$$ container(12): first crop line
% $$$ container(13): last  crop line
% $$$ container(14): first crop pixel
% $$$ container(15): last  crop pixel
% $$$ container(16): mid.crop line
% $$$ container(17): DC_constant  [Hz]
% $$$ container(18): DC_linear    [Hz/s]
% $$$ container(19): DC_quadratic [Hz/s/s]
% $$$ container(20): first crop rsmp line
% $$$ container(21): last  crop rsmp line
% $$$ container(22): first crop rsmp pixel
% $$$ container(23): last  crop rsmp pixel
% $$$ --------------------------------------------------------
% $$$ containerChar(1): crop_file_name
% $$$ containerChar(2): crop_file_format
% $$$ --------------------------------------------------------

container = zeros(23,1); % initalize
containerChar = cell(2,1);

% acquisition start time
mdata = textread(result_file,'%s','delimiter',':');

%%% ovs.factor: AZIMUTH
[container(8),exitCode] = greper('Multilookfactor_azimuth_direction',mdata);
if exitCode == 1;
    container(8) = 1;
else
    container(8) = 1/container(8);
end;

%%% ovs.factor: RANGE
[container(9),exitCode] = greper('Multilookfactor_range_direction',mdata);
if exitCode == 1;
    container(9) = 1;
else
    container(9) = 1/container(9);
end;

% UTC time in secs for first pixel
ind = strmatch('First_pixel_azimuth_time (UTC)', mdata);
ddhr = char(mdata(ind+1)); % date and hour
hms  = [str2num(ddhr(end-1:end)) str2num(char(mdata(ind+2))) ...
        str2num(char(mdata(ind+3)))]; % hh mm ss

container(1) = hms(1)*3600 + hms(2)*60 + hms(3);

% 1way range time to 1st pix in [sec]
container(2) = (greper('Range_time_to_first_pixel (2way) (ms)',mdata)/2)/ 1000;

% PRF [Hz] : accounted for ovsmp factor
[container(3),exitCode] = greper('Pulse_Repetition_Frequency (computed, Hz)',mdata);
if exitCode == 1;
    container(3) = greper('Pulse_Repetition_Frequency (actual, Hz)', mdata);
end
container(3) = container(3) * container(8);

%%% RSR in [Hz]: computed from 2 way time (see definitions) and accounted
%%%              for ovsmp factor
[container(4),exitFlag] = greper('Range_sampling_rate (computed, MHz)',mdata);
if exitFlag == 1;
    container(4) = greper('Range_sampling_rate (leaderfile, MHz)',mdata);
end
container(4) = container(4)*2*container(9)*10^6; % ONE! way in HZ

% $$$ sensor wavelength: dumped in readdata
container(5) = greper('Radar_wavelength',mdata);

%%% lat.center
container(6) = greper('Scene_centre_latitude',mdata);

%%% lon.center
container(7) = greper('Scene_centre_longitude',mdata);

% ---------------------------------------------------

% $$$ number.of.lines[original]: dumped in readdata and accounted for ovsmp
% $$$ factor
container(10)  = greper('Number_of_lines_original',mdata) * container(8);

% $$$ number.of.pixels[original]: dumped in readdata and accounted for ovsmp
% $$$ factor
container(11) = greper('Number_of_pixels_original',mdata) * container(9);

% ---------------------------------------------------

%%% first.crop.line: accounted for ovsmp factor and 'counting'
[container(12),exitCode] = greper('First_line (w.r.t. original_image)',mdata);
if exitCode == 1;
     container(12) = greper('First_line (w.r.t. original_master)',mdata);
end

if container(12) > 1;
    container(12) = container(12)*container(8) - 1; 
end

%%% last.crop.line
[container(13),exitCode] = greper('Last_line (w.r.t. original_image)',mdata);
if exitCode == 1;
    container(13) = greper('Last_line (w.r.t. original_master)',mdata);
end
container(13) = container(13) * container(8);

%%% mid.crop.line
container(16) = (container(12) + container(13))/2;

% ---------------------------------------------------

%%% first.crop.pixel
[container(14),exitCode] = greper('First_pixel (w.r.t. original_image)',mdata);
if exitCode == 1;
    container(14) = greper('First_pixel (w.r.t. original_master)',mdata);
end

if container(14) > 1;
    container(14) = container(14)*container(9) - 1; 
end

%%% last.crop.pixel
[container(15),exitCode] = greper('Last_pixel (w.r.t. original_image)',mdata);
if exitCode == 1;
    container(15) = greper(['Last_pixel (w.r.t. ' 'original_master)'],mdata);
end
container(15)   = container(15) * container(8);

% $$$ Doppler info ----------------------------------------------------

container(17) = greper('Xtrack_f_DC_constant (Hz, early edge)',mdata);
container(18) = greper('Xtrack_f_DC_linear (Hz/s, early edge)',mdata);
container(19) = greper('Xtrack_f_DC_quadratic (Hz/s/s, early edge)',mdata);

% ---------------------------------------------------

if ~isempty(strmatch('*_Start_resample',mdata)) 
  %only when resampling is done

  %%% first.crop.rsmp.line: accounted for ovsmp factor
  [container(20),exitCode] = greper('First_line (w.r.t. original_master)',mdata);
  %if exitCode == 1;
  %     container(20) = greper('First_line (w.r.t. original_master)',mdata);
  %end
  
  %%% last.crop.rsmp.line
  [container(21),exitCode] = greper('Last_line (w.r.t. original_master)',mdata);
  %if exitCode == 1;
  %    container(21) = greper('Last_line (w.r.t. original_master)',mdata);
  %end
  
  
  % ---------------------------------------------------
  
  %%% first.crop.rsmp.pixel
  [container(22),exitCode] = greper('First_pixel (w.r.t. original_master)',mdata);
  %if exitCode == 1;
  %    container(22) = greper('First_pixel (w.r.t. original_master)',mdata);
  %end
  
  %%% last.crop.pixel
  [container(23),exitCode] = greper('Last_pixel (w.r.t. original_master)',mdata);
  %if exitCode == 1;
  %    container(23) = greper(['Last_pixel (w.r.t. ' 'original_master)'],mdata);
  %end
end

% $$$ Crop file info ----------------------------------------------------

% Data_output_file:
containerChar{1} = greper('Data_output_file',mdata,1);

% Data_output_format:
containerChar{2} = greper('Data_output_format',mdata,1);

% $$$ Prec.ORBIT from RAS file ----------------------------------------

ind      = strmatch('NUMBER_OF_DATAPOINTS',mdata);
norb_pts = str2num(char(mdata(ind+1)));
orbit    = str2num(char(mdata(ind+2:ind+1+norb_pts)));

% -----------------------------------------------------
function [value,exitCode] = greper(regexp,stringStream,nonumFlag)
    
    if nargin < 3 || nonumFlag ~= 1;
        nonumFlag='';
        value = 0;
    elseif nonumFlag == 1;
        value = 'null';
    end

    ind = strmatch(regexp,stringStream);
    
    if isempty(ind);

% $$$         error('EXIT: parameter [%s] not in DORIS .result file',regexp);
        exitCode = 1;
    elseif isempty(nonumFlag);
        value = str2num(char(stringStream(ind+1)));
        exitCode = 0;
    else
        value = char(stringStream(ind+1));
        exitCode = 0;
    end
