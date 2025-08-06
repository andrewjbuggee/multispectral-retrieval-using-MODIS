% Collect all neccessary geolocation data for calculating cloud retrievals
% from MODIS

function [sensor,solar,geo] = readMODIS_geolocation(fileName)


    info = hdfinfo(fileName);
    
    

    % load the geolocation data from the MOD03 geolocation hdf file
    % these are the lat-long positions of MODIS pixels on Earths surface
    geo.lat = hdfread(fileName,'Latitude');
    geo.long = hdfread(fileName,'Longitude');

    % Sometimes MODIS won't properly register data. If this happens, it
    % will assign a value of -999 to the latitude and longitude. Let's
    % convert these to nans'

    geo.lat(geo.lat<-90) = nan;
    geo.long(geo.long<-180) = nan;

    % --------------------------------------------------------------------
    % -----------*******------------ NOTE -----------*********------------
    % --------------------------------------------------------------------
    % The azimuth angles are relative to local geodetic north. This is how
    % 0 degrees is defined. Azimuth values range from [-180, 180]. If the
    % measurement is made in the northern hemisphere, the. the positive
    % azimuth values range from due north (0 deg) to due East (90 deg) to
    % due south (180 deg). And moving from due north towards the west would
    % result in negative values from 0 to -180. The same goes for the
    % southern hemisphere. Local geodetic north defines the 0 degree
    % azimuth location. Negative azimuth values are towards local west, and
    % positive azimuth values are towards local east.

    % Read the scale factors for the solar geometry
    solarZenith_scale = info.Vgroup.Vgroup(2).SDS(8).Attributes(4).Value;
    solarAzimuth_scale = info.Vgroup.Vgroup(2).SDS(9).Attributes(4).Value;


    % load solar position data
    solar.azimuth = hdfread(fileName,'SolarAzimuth')*solarZenith_scale; % scale factor included
    solar.zenith = hdfread(fileName,'SolarZenith')*solarAzimuth_scale; % scale factor included
    
    
    % Read the scale factors for the sensor geometry
    sensorZenith_scale = info.Vgroup.Vgroup(2).SDS(5).Attributes(4).Value;
    sensorAzimuth_scale = info.Vgroup.Vgroup(2).SDS(6).Attributes(4).Value;
    sensorRange_scale = info.Vgroup.Vgroup(2).SDS(7).Attributes(4).Value;
    
    % load satellite position data
    % height of the ground location above the Earth ellipsoid (meters)
    sensor.gound_height = double(hdfread(fileName,'Height'));
    % path length from the pixel on the ground to the satellite (meters)
    sensor.range = double(hdfread(fileName,'Range'))*sensorRange_scale; % scale factor included
    sensor.azimuth = double(hdfread(fileName,'SensorAzimuth'))*sensorAzimuth_scale; % scale factor included
    sensor.zenith = double(hdfread(fileName,'SensorZenith'))*sensorZenith_scale; % scale factor included


    % ---- MODIS Parameters -----


end

