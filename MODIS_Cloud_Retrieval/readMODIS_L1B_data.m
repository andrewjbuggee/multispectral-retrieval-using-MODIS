%% ----- READ IN MODIS DATA -----

% this function will read in L1B MODIS data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf


% ---- Description of the different fields within the HDF File -----

%   (1) radiance - units of W/m^2/micron/sr
%   (2) reflectance - units of 1/sr



% By Andrew J. Buggee
%%

function [EV] = readMODIS_L1B_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retrieve the info hdf info structure
info = hdfinfo(fileName);



%% --- Read in a data set at a specifc resolution ---

% First, let's grab two global attributes
% These can be used, along with reflectance, to derive the most precise
% version of radiance (see page 33 of the MODIS L1B user guide)
% Let's read in the Earth-Sun distance as a ratio of 1 AU at the middle
% time of the data granule
%earth_sun_distance = hdfread(fileName, 'Earth-Sun Distance');

% Let's read in the Solar irradiances at 1 AU divided by pi, weighted by a
% detectors relative spectral response for EACH reflective detector
% For each of the 15 1km resolution bands, there are 10 detectors
% For each of the 5 500m resolution bands, there are 20 detectors
% For each of the 2 250m resolution bands, there are 40 detectors
%weighted_solar_irradiance = hdfread(fileName, 'Solar Irradiance on RSB Detectors over pi');





% ----------------------------------------------------------
% Check to see what product we are reading in
% ----------------------------------------------------------

if strcmp(fileName(6:8),'QKM')
    % then we are reading in calibrated Earth View data at 250m resolution
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    uncertainty = 'EV_250_RefSB_Uncert_Indexes';
    
    earthView250 = double(hdfread(fileName,'EV_250_RefSB'));
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    
    EV.radiance = scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m);
    EV.reflectance = scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m);
    
    
    
    
elseif strcmp(fileName(6:8),'HKM')
    
    % --------------------------------------------------------------
    % --------- FOR THIS RETRIEVAL WE ONLY NEED BANDS 1-7 ----------
    % --------------------------------------------------------------
    
    % retreive the scales and offsets for radiance and reflectance. BE
    % CAREFUL!! Each band has a unique scaling and offset. Make sure you have
    % applied the correct scaling. You can check by looking at the info
    % structure
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Radiance scales for bands 3-7 ---
    radianceScales_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(6).Value; % output should be a vector with 5 entries for the bands 3-7
    radianceOffsets_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(7).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % --- Reflectance scales for bands 3-7 ---
    reflectanceScales_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value; % output should be a vector with 5 entries for the bands 3-7
    reflectanceOffsets_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % then we are reading in calibrated Earth View data at 500m resolution,
    % including the 250m resolution bands aggregated to 500m resolution
    
    
    uncertainty = 'EV_500_RefSB_Uncert_Indexes';
    
    earthView250 = double(hdfread(fileName,'EV_250_Aggr500_RefSB'));
    earthView500 = double(hdfread(fileName,'EV_500_RefSB'));
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);
    
    
    
    EV.radiance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m),...
        scalesOffsets2Matrix(m500_scaledIntegers,radianceScales_500m,radianceOffsets_500m));
    
    EV.reflectance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m), ...
        scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales_500m,reflectanceOffsets_500m));
    
    
    
elseif strcmp(fileName(6:8),'1KM')
    % Bands are converted to have 1km resolution
    
    % --------------------------------------------------------------
    % --------- FOR THIS RETRIEVAL WE ONLY NEED BANDS 1-7 ----------
    % --------------------------------------------------------------
    
    % retreive the scales and offsets for radiance and reflectance. BE
    % CAREFUL!! Each band has a unique scaling and offset. Make sure you have
    % applied the correct scaling. You can check by looking at the info
    % structure
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands

    % --- Radiance Uncertainty scales for the first two bands ---
    radianceUncertaintyScales_250m = info.Vgroup.Vgroup(2).SDS(6).Attributes(5).Value; % specfied uncertainty for relfectance in the first two bands
    radianceUncertaintyOffset_250m = info.Vgroup.Vgroup(2).SDS(6).Attributes(6).Value; % specfied uncertainty for relfectance in the first two bands

    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance Uncertainty Scales and Offsets for the first two bands ---
    reflectanceUncertainty_specified_250m = info.Vgroup.Vgroup(2).SDS(6).Attributes(5).Value; % specfied uncertainty for relfectance in the first two bands
    reflectanceUncertainty_scale_250m = info.Vgroup.Vgroup(2).SDS(6).Attributes(6).Value; % specfied uncertainty for relfectance in the first two bands

        
    % --- Radiance scales for bands 3-7 ---
    radianceScales_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(6).Value; % output should be a vector with 5 entries for the bands 3-7
    radianceOffsets_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(7).Value; % output should be a vector with 5 entries for the bands 3-7


    % --- Radiance Uncertainty Scales and Offsets for bands 3-7 ---
    radianceUncertainty_specified_500m = info.Vgroup.Vgroup(2).SDS(9).Attributes(5).Value; % specfied uncertainty for relfectance in the first two bands
    radianceUncertainty_scale_500m = info.Vgroup.Vgroup(2).SDS(9).Attributes(6).Value; % specfied uncertainty for relfectance in the first two bands
    
    % --- Reflectance scales for bands 3-7 ---
    reflectanceScales_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value; % output should be a vector with 5 entries for the bands 3-7
    reflectanceOffsets_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(10).Value; % output should be a vector with 5 entries for the bands 3-7
    
    
    % --- Reflectance Uncertainty Scales and Offsets for bands 3-7 ---
    reflectanceUncertainty_specified_500m = info.Vgroup.Vgroup(2).SDS(9).Attributes(5).Value; % specfied uncertainty for relfectance in the first two bands
    reflectanceUncertainty_scale_500m = info.Vgroup.Vgroup(2).SDS(9).Attributes(6).Value; % specfied uncertainty for relfectance in the first two bands

    
    
    % --------------------------------------------------------------
    % --- NOT USING BANDS 8-36 SO WE DONT NEED EV_1KM_RefSB data ---
    % --------------------------------------------------------------
    
    
    % then we are reading in calibrated Earth View data at 1km resolution
    % including the 250m and 500m data resolution data aggregated to 1km
    % resolution
    
    % Uncertainty Index for 250m EV reflectance product aggregated to 1 km
    % This is the uncertainty for bands 1-2 measured in percent
    uncertainty_250 = hdfread(fileName,'EV_250_Aggr1km_RefSB_Uncert_Indexes');
    
    % 250m Earth Viewing (EV) Science data aggregated to 1 km. This is the
    % data covering the first two MODIS bands
    earthView250 = double(hdfread(fileName,'EV_250_Aggr1km_RefSB')); % first two modis bands aggregated to 1km resolution
    
    
    % Uncertainty Index for 250m EV reflectance product aggregated to 1 km
    % This is the uncertainty for bands 3-7 measured in percent
    uncertainty_500 = hdfread(fileName,'EV_500_Aggr1km_RefSB_Uncert_Indexes');
    
    
    % 500m Earth Viewing (EV) Science data aggregated to 1 km. This is the
    % data covering the MODIS bands 3-7
    earthView500 = double(hdfread(fileName,'EV_500_Aggr1km_RefSB')); % modis bands 3-7 aggregated to 1km resolution
    
    %earthView1000 = hdfread(fileName,'EV_1KM_RefSB'); % modis bands 8 - 36 at 1km resolution
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);

    uncertainty250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(uncertainty_250);   
    uncertainty500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(uncertainty_500);

    EV.radiance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m),...
        scalesOffsets2Matrix(m500_scaledIntegers,radianceScales_500m,radianceOffsets_500m));
    
    EV.reflectance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m), ...
        scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales_500m,reflectanceOffsets_500m));
    
    
    % --- Compute the uncertainty for each band ---
    % This uncertainty is measured as a percent
    % Meaning, its values are between 0 and 100
    EV.reflectanceUncert = cat(3,repmat(reshape(reflectanceUncertainty_specified_250m,1,1,[]),...
        size(uncertainty250_scaledIntegers,1),size(uncertainty250_scaledIntegers,2)).* ...
        exp(uncertainty250_scaledIntegers./repmat(reshape(reflectanceUncertainty_scale_250m,1,1,[]),...
        size(uncertainty250_scaledIntegers,1),size(uncertainty250_scaledIntegers,2))),...
        repmat(reshape(reflectanceUncertainty_specified_500m,1,1,[]),size(uncertainty500_scaledIntegers,1),...
        size(uncertainty500_scaledIntegers,2)).* exp(uncertainty500_scaledIntegers./...
        repmat(reshape(reflectanceUncertainty_scale_500m,1,1,[]),size(uncertainty500_scaledIntegers,1),...
        size(uncertainty500_scaledIntegers,2))));
    
    % --- DONT NEED BANDS 8-36 FOR NOW ---
    
    %     % retrieve the 1km bands, bands 8-36
    %     EV.km1.radiance = scalesOffsets2Matrix(m1000_scaledIntegers,radianceScales,radianceOffsets);
    %     EV.km1.reflectance = scalesOffsets2Matrix(m1000_scaledIntegers,reflectanceScales,reflectanceOffsets);
    %
    %     % the 1 kilometer resolution reflective bands span 8-16 and 26
    %     EV.km1.bands = readMODISbands([8:16,26]);
    
    
elseif strcmp(fileName(6:8),'OBC')
    % then we are reading in the On-Board Calibrator and Engineering Data
    
    
    
    
    
else
    
    error('Filename is not valid! Check to see it is a L1B MODIS data file')
    
end




% ----------------------------------------------------------
% ------------ Read in MODIS L1B Swath Meta Data -----------
% ----------------------------------------------------------

L1B_metadata = hdfread(fileName, 'Level 1B Swath Metadata');

raw_time = L1B_metadata{5};     % TAI time (number of sec. since 1/1/93)


% there are typically 203 scans for each granule. Since there are 10
% detectors for the 1KM data, that means there are 2030 total measurements
% in the along track direction. Along the scan direction, there are 1354
% measurements, or pixels. 
% Let's interpolate the get the time of each pixel in the along-scan
% direction

modis_pixel_time_along_scan = zeros(size(EV.radiance, 2), length(raw_time));

for nn=1:(length(raw_time)-1)
    
    time_from_one_scan_to_the_next = linspace(raw_time(nn), raw_time(nn+1), size(EV.radiance,2)+1);
    modis_pixel_time_along_scan(:,nn) = time_from_one_scan_to_the_next(1:end-1)';

end

% the last column is the final scan time
modis_pixel_time_along_scan(:,end) = linspace(raw_time(end), raw_time(end) + (time_from_one_scan_to_the_next(end-1) - ...
    time_from_one_scan_to_the_next(1)), size(EV.radiance, 2))';

% to get the time of every pixel in a full 1KM resolution granule, we
% repeat values in the along-track dimension. For 1KM resolution data,
% there are 10 sensors. So every column (each scan) should be repeated 10
% times.
modis_pixel_time = [];

for nn = 1:length(raw_time)
    
    modis_pixel_time = [modis_pixel_time, repmat(modis_pixel_time_along_scan(:,nn), 1, 10)];

end

% now let's convert these times to UTC
% MODIS EV Sector time starts on 1 Jan 1993
epoch = datetime(1993,1,1,'TimeZone','UTCLeapSeconds');

% now lets add secods according to the EV sector start time measurements
EV.pixel_time_UTC = (epoch + seconds(modis_pixel_time))';

% let's also store the time in decimal foramt, where the decimal represents
% the fraction of the hour (e.g. 14.5 is 14:30 UTC)
EV.pixel_time_decimal = EV.pixel_time_UTC.Hour + EV.pixel_time_UTC.Minute/60 + EV.pixel_time_UTC.Second/3600;





end





