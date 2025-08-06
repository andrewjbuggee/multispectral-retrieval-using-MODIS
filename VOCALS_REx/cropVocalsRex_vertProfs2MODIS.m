

function vocalsRex = cropVocalsRex_vertProfs2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold, modis, modisInputs)


% ----- Find all vertical profiles within VOCALS-REx data ------
%vert_prof = find_verticalProfiles_VOCALS_REx(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold);
vert_prof = find_verticalProfiles_VOCALS_REx_ver2(vocalsRex, lwc_threshold, stop_at_max_lwc,...
                        Nc_threshold, modisInputs.which_computer);


%% Let's step through each vertical profile and find the MODIS pixel that overlaps


% --------------------------------------------------------------------
% ------ Find the MODIS pixels to use for the vertical retrieval -----
% --------------------------------------------------------------------

% store the MODIS latitude and longitude
modis_lat = modis.geo.lat;
modis_long = modis.geo.long;

% store modis pixel time
modis_pixel_time = modis.EV1km.pixel_time_decimal;


% we will be computing the arclength between points on an ellipsoid
% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with units of meters.
wgs84 = wgs84Ellipsoid("m");

% Set up an empty array for each vocals-rex profile
modis_minDist = zeros(1, length(vert_prof));
time_diff_MODIS_VR = zeros(1, length(vert_prof));
within_5min_window = zeros(1, length(vert_prof));

% Step through each vertical profile and find the MODIS pixel that overlaps
% with the mid point of the in-situ sample
parfor nn = 1:length(vert_prof)


    dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, vert_prof(nn).latitude(round(end/2)),...
        vert_prof(nn).longitude(round(end/2)), wgs84);

    [modis_minDist(nn), index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');            % m - minimum distance


    % compute the time between the modis pixel closest to the VR sampled
    % path and the time VOCALS was recorded
    time_diff_MODIS_VR(nn) = abs(modis_pixel_time(index_minDist) - ...
        vert_prof(nn).time_utc(round(end/2))) * 60;                         % minutes

    % A single MODIS granule takes 5 minutes to be collected. First, lets
    % determine if any profiles were recorded within the the 5 minute window
    % when MODIS collected data.

    within_5min_window(nn) = any(vert_prof(nn).time_utc >= modis_pixel_time(1,1) & ...
        vert_prof(nn).time_utc <= modis_pixel_time(end,end));

end

%% Check the time difference and the distance between the closest pixel and VR


% First, keep all profiles recorded within the MODIS granule time window
% and have a minimum distance with the closest MODIS pixel of less than 1km
if any(within_5min_window==true & modis_minDist<=1000)
    
    index_vertProfs_2keep = find(within_5min_window==true & modis_minDist<=1000);

% Next, keep profiles with the smallest time difference as long as the
% distance between the closest MODIS pixel and the VR profile is less than
% 1km
elseif any(time_diff_MODIS_VR<5 & modis_minDist<1000)
    
%     % First, check to see if this is the MODIS file from 11 Nov. 2008 recorded at 1850 UTC.
%     % This file is tricky because two vertical profiles recorded right
%     % before and right after the MODIS granule are very different. The
%     % first has an adiabatic droplet profile, the second is nearly
%     % homogeneous.
%     if strcmp(modisInputs.L1B_filename, 'MYD021KM.A2008316.1850.061.2018039033053.hdf')==true
%         % pick the second closest profile
%         [~, idx_min] = min(time_to_halfwayPoint);
%         % get rid of this value
%         time_to_halfwayPoint(idx_min) = inf;
%         % grab the next vertical profile closest in time
%         [~, index_vertProfs_2keep] = min(time_to_halfwayPoint);
%     else
% 
%         [~, index_vertProfs_2keep] = min(time_to_halfwayPoint);
% 
%     end


    % save vertical profiles that were recorded within 5 minutes of the
    % MODIS recording and the distance between VR and the closest pixel is
    % less than 1km
    index_vertProfs_2keep = find(time_diff_MODIS_VR<5 & modis_minDist<1000);



else

    % If there isn't a profile sampled within the 5 minute window, or
    % within 5 minutes of a MODIS pixel, take the profile with the shortest
    % amount of time between it's recording and MODIS, as long as it's
    % within 1 km of the pixel in question
    [~, min_time_idx] = min(time_diff_MODIS_VR);
    if modis_minDist(min_time_idx)<1000

        index_vertProfs_2keep = min_time_idx;

    else

        error([newline, 'The minimum time difference between a VR in-situ profile and a MODIS pixel ',...
            'was ', num2str(min(time_diff_MODIS_VR)), ' minutes, and the distance to the closest pixel ',...
            'was ', num2str(modis_minDist(min_time_idx)), ' meters.', newline])
    end


end

% I dont' know what to do if there is more than 1 profile that satisfies
% the above criteria
if length(index_vertProfs_2keep)>1
    error([newline, "Code isn't set up for more than 1 profile.", newline])
end

% Let's keep the vertical profile that is closest to MODIS
clear vocalsRex;

% Keep only the profiles listed in index_vertProfs_2keep
vocalsRex = vert_prof(index_vertProfs_2keep);



%%
% --------------------------------------------------------------------
% ------ Find the MODIS pixels to use for the vertical retrieval -----
% --------------------------------------------------------------------


% Set up an empty array for each vocals-rex data point
modisIndex_minDist = zeros(1, length(vocalsRex.latitude));
modis_minDist = zeros(1, length(vocalsRex.latitude));

% Store Vocals-Rex lat and long
vr_lat = vocalsRex.latitude;
vr_long = vocalsRex.longitude;

% store length of data points for VOCALS-REx
n_data_VR = length(vr_lat);


% First, let's find the MODIS pixel closest to ALL VOCALS-REx locations
% throughout the sampled vertical profile
parfor nn = 1:n_data_VR


    dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, vr_lat(nn), vr_long(nn), wgs84);

    [modis_minDist(nn), modisIndex_minDist(nn)] = min(dist_btwn_MODIS_and_VR, [], 'all');       % m - minimum distance between VR and closest MODIS pixel


end

% Did the VOCALS-REx measurement take place before or after the MODIS
% sample? Use the pixel closest to VOCALS-REx
[global_min_dist, min_idx] = min(modis_minDist);
VR_sampled_before_MODIS = vocalsRex.time_utc(min_idx) < modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx));


if modisInputs.flags.useAdvection==false

    % ------------------------ NO ADVECTION -------------------------
    % if advection flag is false, simply find the MODIS pixels closest in
    % space to the vocals-rex in-situ measurements




else


    % ----------------------- USE ADVECTION -----------------------
    % if advection is true, use the measured windspeed and direction to
    % project where the cloud was at the time of the MODIS overpass



    % project EACH vocalsRex data point using the measured wind speed and
    % direction
    horz_wind_speed = reshape(vocalsRex.horz_wind_speed,[], 1);    % m/s

    % use EACH wind direction
    % This is the direction the wind is coming from, NOT the direction the
    % wind is blowing torwards
    wind_from_direction = vocalsRex.horz_wind_direction;         % degrees from north (0 deg)

    % compute the direction the wind is blowing towards using modulo
    % arithmetic
    wind_direction = mod(wind_from_direction + 180, 360);               % degrees from north (0 deg)

    % If VOCALS-REx sampled after the MODIS pixel closest to the in-situ
    % measured path was recorded, then in order to line up the same cloudy
    % region for comparison, we need to move the VOCALS-REx flight backwards
    % in time along the direction of where the wind is comming from. That's
    % because the we want to use the cloud that moved into the VOCALS-REx
    % position once VOCALS-REx made it's measurement

    % compute the time in seconds between the MODIS overpass and the
    % vocalsRex in-situ measurement
    d_time_sec = abs(vocalsRex.time_utc(min_idx) - modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx)))*3600;           % sec


    % compute the horizontal distance travelled by the cloud during this
    % time
    dist_m = horz_wind_speed.*d_time_sec;                % meters travelled


    % First, find whether or not MODIS passed overhead before or after the
    % VOCALS-REX made in-situ measurements.
    % if true, then MODIS passed overhead before VOCALS sampled the cloud
    % set this to be a value of 1
    modis_before_vocals = modis.EV1km.pixel_time_decimal(modisIndex_minDist(min_idx))<vocalsRex.time_utc;


    % if position change is a logical 1, we have to move the VOCALS-REx
    % in-situ cloud position backward in time into the direction the
    % wind is coming from
    azimuth_angle = zeros(n_data_VR,1);
    azimuth_angle(modis_before_vocals) = wind_from_direction(modis_before_vocals);

    % if position change is logical 0, we have to move the VOCALS-REx
    % in-situ cloud position forwards in time. This time we move the data
    % along the direction the wind is moving into
    azimuth_angle(~modis_before_vocals) = wind_direction(~modis_before_vocals);


    % Use the reckon function to compute the new lat long position.
    % Inputs are the original lat/long, the distance travelled and the
    % azimuth, which is the angle between the direction of travel and
    % true north. Or, as MATLAB puts it, it is the angle between the
    % local meridian line and the direction of travel, where the angle
    % is swept out in the clockwise direction.
    % (0 - due north, 90 - due east, 180 - due west etc.)

    [lat_withAdvection, long_withAdvection] = reckon(vr_lat', vr_long', dist_m, azimuth_angle, wgs84);

    
    % Set up an empty array for each vocals-rex data point
    modisIndex_minDist = zeros(1, length(vocalsRex.latitude));
    modis_minDist = zeros(1, length(vocalsRex.latitude));
    parfor nn = 1:n_data_VR

        dist_btwn_MODIS_and_VR = distance(modis_lat, modis_long, lat_withAdvection(nn), long_withAdvection(nn), wgs84);

        [min_dist, index_minDist] = min(dist_btwn_MODIS_and_VR, [], 'all');
        % save this index so we know what MODIs pixel to use in comparison
        modisIndex_minDist(nn) = index_minDist;
        modis_minDist(nn) = min_dist;            % meters

    end

    % Save the new lat and long
    vocalsRex.lat_withAdvection = lat_withAdvection;
    vocalsRex.long_withAdvection = long_withAdvection;




end




% Store the modis index values and the distance from VOCALS to the pixel
vocalsRex.modisIndex_minDist = modisIndex_minDist;
vocalsRex.modis_minDist = modis_minDist;            % meters




% ---- Check to see if all indices are unique. If not, delete the
% redundancy

[~, idx_unique] = unique(vocalsRex.modisIndex_minDist);
idx_unique_logic = ismember(1:length(vocalsRex.modisIndex_minDist), idx_unique);

% delete the indices that are not found above

vocalsRex.modisIndex_minDist = vocalsRex.modisIndex_minDist(idx_unique_logic);

% compute the time difference in minutes between the VOCALS_REx profile and
% each MODIS pixel selected
vocalsRex.timeDiff_withMODIS = abs(modis_pixel_time(vocalsRex.modisIndex_minDist) - ...
        vocalsRex.time_utc(round(end/2))) * 60;                         % minutes








end