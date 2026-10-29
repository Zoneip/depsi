function [container, orbit] = metadata(result_file)

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
% $$$ [23.04.08, 14:36]: different CONTAINER vector structure

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
% $$$ --------------------------------------------------------

container = zeros(1,19); % initalize

% acquisition start time
mdata = textread(result_file,'%s','delimiter',':');

%%% ovs.factor: AZIMUTH
ind = strmatch('Multilookfactor_azimuth_direction',mdata);
if isempty(ind);
    container(8) = 1;
else
    container(8) = 1/str2num(char(mdata(ind+1)));
end;

%%% ovs.factor: RANGE
ind = strmatch('Multilookfactor_range_direction',mdata);
if isempty(ind);
    container(9)=1;
else
    container(9) = 1/str2num(char(mdata(ind+1)));
end;

% UTC time in secs for first pixel
ind = strmatch('First_pixel_azimuth_time (UTC)', mdata);
ddhr = char(mdata(ind+1)); % date and hour

hms  = [str2num(ddhr(end-1:end)) str2num(char(mdata(ind+2))) ...
        str2num(char(mdata(ind+3)))]; % hh mm ss

container(1) = hms(1)*3600 + hms(2)*60 + hms(3);

% 1way range time to 1st pix in [sec]
ind      = strmatch('Range_time_to_first_pixel (2way) (ms)',mdata);
container(2) = (str2num(char(mdata(ind+1)))/2)/1000; % one way in seconds

% PRF [Hz]
ind      = strmatch('Pulse_Repetition_Frequency (computed, Hz)',mdata);
if isempty(ind);
    ind      = strmatch('Pulse_Repetition_Frequency (actual, Hz)', ...
                        mdata);
end
container(3) = str2num(char(mdata(ind+1))) * container(8);

%%% RSR in [Hz]: computed from 2 way time, see definitions
ind      = strmatch('Range_sampling_rate (computed, MHz)',mdata);
if isempty(ind);
    ind      = strmatch('Range_sampling_rate (leaderfile, MHz)',mdata);
end
container(4) = (str2num(char(mdata(ind+1)))*1e6)*2*container(9); % ONE! way in HZ

% $$$ sensor wavelength: dumped in readdata
ind = strmatch('Radar_wavelength',mdata); %ind = ind(end);
container(5) = str2num(char(mdata(ind+1))); 

%%% lat.center
ind      = strmatch('Scene_centre_latitude',mdata);
container(6) = str2num(char(mdata(ind+1)));

%%% lon.center
ind      = strmatch('Scene_centre_longitude',mdata);
container(7) = str2num(char(mdata(ind+1)));

% ---------------------------------------------------

% $$$ number.of.lines[original]: dumped in readdata
ind = strmatch('Number_of_lines_original',mdata); %ind = ind(end);
container(10) = str2num(char(mdata(ind+1))) * container(8); 

% $$$ number.of.pixels[original]: dumped in readdata
ind = strmatch('Number_of_pixels_original',mdata); %ind = ind(end);
container(11) = str2num(char(mdata(ind+1))) * container(9);

% ---------------------------------------------------

%%% first.crop.line
ind      = strmatch('First_line (w.r.t. original_image)',mdata);
if isempty(ind);
    ind      = strmatch('First_line (w.r.t. original_master)',mdata);
end

first_l  = str2num(char(mdata(ind+1))) * container(8);

if container(7) > 1;
    container(12) = first_l - 1; 
else
    container(12) = first_l; 
end

%%% last.crop.line
ind      = strmatch('Last_line (w.r.t. original_image)',mdata);
if isempty(ind);
    ind      = strmatch('Last_line (w.r.t. original_master)',mdata);
end
last_l   = str2num(char(mdata(ind+1)))  * container(8);
container(13) = last_l;

%%% mid.crop.line
% $$$ container(9) = round((first_l + last_l)/2);
container(16) = (first_l + last_l)/2;

% ---------------------------------------------------

%%% first.crop.pixel
ind      = strmatch('First_pixel (w.r.t. original_image)',mdata);
if isempty(ind);
    ind      = strmatch('First_pixel (w.r.t. original_master)',mdata);
end
first_px = str2num(char(mdata(ind+1))) * container(9);

if container(8) > 1;
    container(14) = first_px - 1; 
else
    container(14) = first_px; 
end

%%% last.crop.pixel
ind      = strmatch('Last_pixel (w.r.t. original_image)',mdata);
if isempty(ind);
    ind      = strmatch('Last_pixel (w.r.t. original_master)',mdata);
end
last_px   = str2num(char(mdata(ind+1)))  * container(8);
container(15) = last_px;


% $$$ Doppler info ----------------------------------------------------

ind = strmatch('Xtrack_f_DC_constant (Hz, early edge)',mdata);
container(17) = str2num(char(mdata(ind+1)));

ind = strmatch('Xtrack_f_DC_linear (Hz/s, early edge)',mdata);
container(18) = str2num(char(mdata(ind+1)));

ind = strmatch('Xtrack_f_DC_quadratic (Hz/s/s, early edge)',mdata);
containery(19) = str2num(char(mdata(ind+1)));


% $$$ Prec.ORBIT from RAS file ----------------------------------------

ind      = strmatch('NUMBER_OF_DATAPOINTS',mdata);
norb_pts = str2num(char(mdata(ind+1)));
orbit    = str2num(char(mdata(ind+2:ind+1+norb_pts)));

% EOF
