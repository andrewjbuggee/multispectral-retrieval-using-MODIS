%% Finding VOCALS-REx data where the plane flies cleanly through a cloud layer
% Save these vertical profiles


% INPUTS:
% -------
%       (1) vocals-rex data set

%       (2) LWC_threshold - (g/m^3) this is a threshold that helps clean
%       the data. If set to 0, the vertical profiles will have values on
%       either end of the true vertical droplet profile that represent
%       measurements outside of what we want. the LWC_threshold can help
%       clean the data by limiting the data to have a certain liquid water
%       content threshold. For example. Painemal and Zuidema (2011) defined
%       the cloud top boundary as the level where the LWC was 0.03 g/m^3.
%       Beyond this, the excluded the data. This value will be used to
%       truncate data before and after the vertical profile.


%       (3) stop_at_max_lwc - this is a logical flag, whose values can either
%       be true or false. If true, each vertical profile found will be
%       truncated at the peak value of the liquid water content. Most
%       vertical droplet profiles will have a LWC that grows with altitude.
%       If this is true, the profile usually peters out shortly after the
%       peak LWC value. IF this is false, the data will be kept until the
%       the LWC is below the LWC threshold defined above.


%       (4) Nc_threshold - (droplets/cm^3) this is a threshold that helps
%       the function find profiles with some tangible physical meaning. At
%       times there are confounding measurements where the LWC is greater
%       than the defined threshold but the total number concentration is
%       less than 1. Typically this coincides with erroneous droplet size
%       estimates. This value will be used to ensure only contiguous data
%       above this threshold will satisfy the profile search.



function [vert_profs] = find_verticalProfiles_VOCALS_REx_ver2(vocalsRex, LWC_threshold, stop_at_max_lwc, Nc_threshold, which_computer)


% Create an new structure for the vertical profiles
% grab fieldnames from VocalsRex
fields = fieldnames(vocalsRex);


% ------------------------------------------------------------
% ------------------ Vertical Profile Requirements -----------
% ------------------------------------------------------------
% dz/dt must be non-zero.
% Total Nc has to start at a value below 1
% Total Nc has to end at a value below 1
% LWC has to start below threshold
% LWC has to end below threshold
% All intervening data points must exceed LWC and Nc thresholds
% -------------------------------------------------------------

% First lets crop the data only to those portions where the plane is
% ascending or descending

dz_dt = diff(vocalsRex.altitude)./diff(vocalsRex.time);            % m/s

% Compute the mean with a sliding window for every n data points
% This will smooth out the data and make the horizontal flight segments
% easier to find. It's important we ignore the horizontal segments and find
% only the vertical profiles

n_window = 20;
% append a zero to the end so it has the same length as the original time
% vector
dz_dt_mean = [movmean(dz_dt,n_window), 0];


% Looking at a plot of time versus dz/dt as well as the plane's altitude,
% it's clear the plane's vertical velocity exceeds +/- 2 m/s on average when
% it ascends or descends, respectively. Use this to seperate the data
% between profiles measured via ascending or descending.
vertical_velocity_threshold = 2;    % m/s


% Find regions where total Nc transitions from less than 1 to greater than
% 1, LWC transitions from below 0.03 to greater than 0.03, and dz/dt~=0
idx_dz_dt = abs(dz_dt_mean)>vertical_velocity_threshold;



% at what idx's does the total Nc change from below 1 to above 1
% these indexes bookend vertical profiles. They should happen before and
% after.
idx_Nc_lwc_transition = [0, diff(vocalsRex.total_Nc > Nc_threshold & vocalsRex.lwc > LWC_threshold)];

% profiles with start with a transition idx of 1 and end with a transition
% idx of -1. Search in between these indexes, where the separation is
% atleast 3 indexes
% first find all transition idxs with a value of 1
idx_1 = find(idx_Nc_lwc_transition==1);

% define minimum number of points to keep as a profile
length_threshold = 15;

% define the number of indexes to check before and after profile's first
% and last index
num_ba = 20;

% If needed use more indices before and after the profile
num_ba_long = 30;


% Start with zero profiles and add a 1 for each profile found
profile_num = 0;



% step through each idx
for nn = 1:length(idx_1)

    % For every cloud idx_1 index where the we transition from an
    % atmosphere with a total number concentration and liquid water content
    % below the defined thresholds to a cloud with values above the
    % thresholds, find the next index where both the total
    % number concentration AND the LWC are below the thresholds
    idx_cloud_boundary = find(vocalsRex.total_Nc(idx_1(nn)+1:end) < Nc_threshold & ...
        vocalsRex.lwc(idx_1(nn)+1:end) < LWC_threshold);

    % the index below is the LAST value before BOTH total number
    % concentration and liquid water content dip below their respective
    % thresholds
    idx_cloud_boundary = idx_cloud_boundary(1) + idx_1(nn) -1;

    % Each time the total number concentration and LWC thresholds are
    % exceed, our idx_Nc_LWC_transition will be marked with a 1. When
    % either of these thresholds are no longer met, our
    % idx_Nc_LWC_transition will be mark a -1. So every single 1 will have
    % a corresponding -1
    % Search between each [1,-1] pair and ensure while the intervening
    % measurements took place, the mean dz/dt of the plane was above some
    % threshold. This indicates that plane was eitehr climbing or
    % descending during the entire cloud that was sampled. Lets also ensure
    % that these profiles have a length of atleast our length_threshold

    % is the plane descending or ascending between the cloud boundaries?
    %     ascend_or_descend_throughout_cloud = all(idx_dz_dt(idx_1(nn):idx_minus1(nn)));
    ascend_or_descend_throughout_cloud = all(idx_dz_dt(idx_1(nn):idx_cloud_boundary));

    % is the measured cloud segment longer than our threshold?
    %     meets_length_requirement = length(idx_1(nn):idx_minus1(nn))>=length_threshold;
    meets_length_requirement = length(idx_1(nn):idx_cloud_boundary)>length_threshold;

    % check to make sure the median value of the 7 data points before and
    % after are below the thresholds. We use the median because the range
    % of values can be in the hundreds over a few data points, skewing the
    % average significantly. For LWC we can use the average since the range
    % is much smaller

    % check if idx_1(nn)-num_ba is a positive number
    if (idx_1(nn)-num_ba)>0

        before_profile_below_thresholds = median(vocalsRex.total_Nc(idx_1(nn)-num_ba:idx_1(nn)-1)) < Nc_threshold & ...
            mean(vocalsRex.lwc(idx_1(nn)-num_ba:idx_1(nn)-1)) < LWC_threshold;

    else
        % Start from index 1
        before_profile_below_thresholds = median(vocalsRex.total_Nc(1:idx_1(nn)-1)) < Nc_threshold & ...
            mean(vocalsRex.lwc(1:idx_1(nn)-1)) < LWC_threshold;
    end

    % check if idx_cloud_boundary+num_ba>length of data set
    if (idx_cloud_boundary+num_ba)<length(vocalsRex.time)

        after_profile_below_thresholds = median(vocalsRex.total_Nc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba)) < Nc_threshold & ...
            mean(vocalsRex.lwc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba)) < LWC_threshold;

    else
        % end at the last index
        after_profile_below_thresholds = median(vocalsRex.total_Nc(idx_cloud_boundary+1:end)) < Nc_threshold & ...
            mean(vocalsRex.lwc(idx_cloud_boundary+1:end)) < LWC_threshold;
    end

    % The above logical statements tend discard layered cloud systems where
    % there is a very short decrease in the total number concentration and
    % liquid water content, typically over just a few data points (seconds)
    % and then the values increase above the thresholds again.


    %idx_1(nn)
    % if both are true, keep profile
    if ascend_or_descend_throughout_cloud==true && meets_length_requirement==true &&...
            before_profile_below_thresholds==true && after_profile_below_thresholds==true

        % grab those indices and store it as a vertical profile
        profile_num = profile_num+1;

        vert_profs(profile_num) = vocalsRex;

        for ff = 1:length(fields)

            % now remove all data points outside of the vertical profile

            if numel(vert_profs(profile_num).(fields{ff}))==length(idx_dz_dt)

                % all the time data will have be a vector with the same
                % length as our calculation of vertical velocity

                % reshape all fields so that time increases with increasing
                % column number (row vector)
                vert_profs(profile_num).(fields{ff}) = reshape(vert_profs(profile_num).(fields{ff})(idx_1(nn):idx_cloud_boundary),...
                    1, []);


            elseif numel(vert_profs(profile_num).(fields{ff}))>length(idx_dz_dt)
                % only one field is a matrix with more values than the time
                % vector, and thats the matrix for the size distribution,
                % with rows representing different size bins and columns
                % are along the time dimension
                vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff})(:, idx_1(nn):idx_cloud_boundary);

            elseif numel(vert_profs(profile_num).(fields{ff}))<length(idx_dz_dt)

                % some miscallaneous fields have non-timed information
                % keep all this information
                vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff});

            end


        end


        % Let's check to see if there is a short segment over which the
        % measured values drop below the defined thersholds, potentially
        % indicated a multilayered cloud
    elseif ascend_or_descend_throughout_cloud==true && meets_length_requirement==true &&...
            before_profile_below_thresholds==true && after_profile_below_thresholds==false


        % In this scenario, there is a cloud with a short segement where
        % the measured values decrease below the defined thresholds

        % Find the next index where the values increase above the defined
        % thersholds
        idx_cloud_boundary2 = find(vocalsRex.total_Nc(idx_cloud_boundary+1:end) > Nc_threshold & ...
            vocalsRex.lwc(idx_cloud_boundary+1:end) > LWC_threshold);

        % the index below is the next value where BOTH total number
        % concentration and liquid water content are greater than their
        % respective thresholds
        idx_cloud_boundary2 = idx_cloud_boundary2(1) + idx_cloud_boundary;

        % compute the number of measurements that drop below the define
        % thresholds
        length_between_multi_layers = idx_cloud_boundary2 - idx_cloud_boundary - 1;

        % Check to see for how many data points the next cloud layer
        % satisfies both thresholds
        idx_cloud_boundary3 = find(vocalsRex.total_Nc(idx_cloud_boundary2+1:end) < Nc_threshold & ...
            vocalsRex.lwc(idx_cloud_boundary2+1:end) < LWC_threshold);

        % the index below is the next value where BOTH total number
        % concentration and liquid water content are greater than their
        % respective thresholds
        idx_cloud_boundary3 = idx_cloud_boundary3(1) + idx_cloud_boundary2 - 1;

        % How many measurements make up the second cloud layer?
        length_of_second_layer = idx_cloud_boundary3 - idx_cloud_boundary2;

        % Compute the median total number concentration and the mean lwc
        % over a longer segment of data after the end of the profile. Start
        % at the end of the first boundary found. Use the median for boths
        % because the set of values can span two orders of magnitude
        after_profile_below_thresholds_long = median(vocalsRex.total_Nc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba_long)) < Nc_threshold & ...
            median(vocalsRex.lwc(idx_cloud_boundary+1:idx_cloud_boundary+num_ba_long)) < LWC_threshold;


        % Make sure both total number concentration and liquid water
        % content are below the defined thresholds after the second cloud
        % layer
        after_profile_below_thresholds2 = median(vocalsRex.total_Nc(idx_cloud_boundary3:idx_cloud_boundary3+num_ba)) < Nc_threshold & ...
            mean(vocalsRex.lwc(idx_cloud_boundary3:idx_cloud_boundary3+num_ba)) < LWC_threshold;


        % If the number of measurements for the second layer is greater
        % than the number of measurements for the gap between the two
        % layers, and the measurements of Nc and LWC are below the
        % thresholds after the second layer for atleast num_ba
        % measurements, accept the entire layer as a cloud.
        if length_between_multi_layers<5 && after_profile_below_thresholds_long==true &&...
                after_profile_below_thresholds2==true

            % If these are all met, then let's accept the full profile
            % grab those indices and store it as a vertical profile
            profile_num = profile_num+1;

            vert_profs(profile_num) = vocalsRex;

            for ff = 1:length(fields)

                % now remove all data points outside of the vertical profile

                if numel(vert_profs(profile_num).(fields{ff}))==length(idx_dz_dt)

                    % all the time data will have be a vector with the same
                    % length as our calculation of vertical velocity

                    % reshape all fields so that time increases with increasing
                    % column number (row vector)
                    vert_profs(profile_num).(fields{ff}) = reshape(vert_profs(profile_num).(fields{ff})(idx_1(nn):idx_cloud_boundary3),...
                        1, []);


                elseif numel(vert_profs(profile_num).(fields{ff}))>length(idx_dz_dt)
                    % only one field is a matrix with more values than the time
                    % vector, and thats the matrix for the size distribution,
                    % with rows representing different size bins and columns
                    % are along the time dimension
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff})(:, idx_1(nn):idx_cloud_boundary3);

                elseif numel(vert_profs(profile_num).(fields{ff}))<length(idx_dz_dt)

                    % some miscallaneous fields have non-timed information
                    % keep all this information
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff});

                end


            end


        end



        % Let's check to see if there is a short segment over which the
        % measured values drop below the defined thersholds, potentially
        % indicated a multilayered cloud
    elseif ascend_or_descend_throughout_cloud==true && meets_length_requirement==true &&...
            before_profile_below_thresholds==false && after_profile_below_thresholds==true


        % In this scenario, there is a cloud with a short segement where
        % the measured values decrease below the defined thresholds

        % Find the previous index where the values increase above the defined
        % thersholds
        idx_cloud_boundary0 = find(vocalsRex.total_Nc(1:idx_1(nn)-1) > Nc_threshold & ...
            vocalsRex.lwc(1:idx_1(nn)-1) > LWC_threshold);

        % the index below is the previous value where BOTH total number
        % concentration and liquid water content are greater than their
        % respective thresholds
        idx_cloud_boundary0 = idx_cloud_boundary0(end);

        % compute the number of measurements that drop below the define
        % thresholds
        length_between_multi_layers = idx_1(nn) - idx_cloud_boundary0;

        % Check to see for how many data points the next cloud layer
        % satisfies both thresholds
        idx_cloud_boundary_minus1 = find(vocalsRex.total_Nc(1:idx_cloud_boundary0-1) < Nc_threshold & ...
            vocalsRex.lwc(1:idx_cloud_boundary0-1) < LWC_threshold);

        % the index below is the next value where BOTH total number
        % concentration and liquid water content are greater than their
        % respective thresholds
        idx_cloud_boundary_minus1 = idx_cloud_boundary_minus1(end);

        % How many measurements make up the second cloud layer?
        length_of_second_layer = idx_cloud_boundary0 - idx_cloud_boundary_minus1;

        % Compute the median total number concentration and the mean lwc
        % over a longer segment of data after the start of the profile. Start
        % at the end of the first boundary found. Use the median for both
        % because the set of values can span two orders of magnitude
        
        before_profile_below_thresholds_long = median(vocalsRex.total_Nc(idx_1(nn)-num_ba_long:idx_1(nn)-1)) < Nc_threshold & ...
            median(vocalsRex.lwc(idx_1(nn)-num_ba_long:idx_1(nn)-1)) < LWC_threshold;

        % Make sure both total number concentration and liquid water
        % content are below the defined thresholds after the second cloud
        % layer
        before_profile_below_thresholds2 = median(vocalsRex.total_Nc(idx_cloud_boundary_minus1-num_ba_long:idx_cloud_boundary_minus1)) < Nc_threshold & ...
            mean(vocalsRex.lwc(idx_cloud_boundary_minus1-num_ba_long:idx_cloud_boundary_minus1)) < LWC_threshold;


        % If the number of measurements for the second layer is greater
        % than the number of measurements for the gap between the two
        % layers, and the measurements of Nc and LWC are below the
        % thresholds after the second layer for atleast num_ba
        % measurements, accept the entire layer as a cloud.
        if length_between_multi_layers<5 && before_profile_below_thresholds_long==true &&...
                before_profile_below_thresholds2==true

            % If these are all met, then let's accept the full profile
            % grab those indices and store it as a vertical profile
            profile_num = profile_num+1;

            vert_profs(profile_num) = vocalsRex;

            for ff = 1:length(fields)

                % now remove all data points outside of the vertical profile

                if numel(vert_profs(profile_num).(fields{ff}))==length(idx_dz_dt)

                    % all the time data will have be a vector with the same
                    % length as our calculation of vertical velocity

                    % reshape all fields so that time increases with increasing
                    % column number (row vector)
                    vert_profs(profile_num).(fields{ff}) = reshape(vert_profs(profile_num).(fields{ff})(idx_cloud_boundary_minus1:idx_cloud_boundary),...
                        1, []);


                elseif numel(vert_profs(profile_num).(fields{ff}))>length(idx_dz_dt)
                    % only one field is a matrix with more values than the time
                    % vector, and thats the matrix for the size distribution,
                    % with rows representing different size bins and columns
                    % are along the time dimension
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff})(:, idx_cloud_boundary_minus1:idx_cloud_boundary);

                elseif numel(vert_profs(profile_num).(fields{ff}))<length(idx_dz_dt)

                    % some miscallaneous fields have non-timed information
                    % keep all this information
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff});

                end


            end


            % If the gap inbetween the two cloud layers is larger than the
            % second layer, let's only accept the original cloud layer
        elseif length_between_multi_layers>length_of_second_layer && before_profile_below_thresholds_long==true &&...
                before_profile_below_thresholds2==true


            % If theses conditions are met, let's keep the original profile
            % indexes
            profile_num = profile_num+1;

            vert_profs(profile_num) = vocalsRex;

            for ff = 1:length(fields)

                % now remove all data points outside of the vertical profile

                if numel(vert_profs(profile_num).(fields{ff}))==length(idx_dz_dt)

                    % all the time data will have be a vector with the same
                    % length as our calculation of vertical velocity

                    % reshape all fields so that time increases with increasing
                    % column number (row vector)
                    vert_profs(profile_num).(fields{ff}) = reshape(vert_profs(profile_num).(fields{ff})(idx_1(nn):idx_cloud_boundary),...
                        1, []);


                elseif numel(vert_profs(profile_num).(fields{ff}))>length(idx_dz_dt)
                    % only one field is a matrix with more values than the time
                    % vector, and thats the matrix for the size distribution,
                    % with rows representing different size bins and columns
                    % are along the time dimension
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff})(:, idx_1(nn):idx_cloud_boundary);

                elseif numel(vert_profs(profile_num).(fields{ff}))<length(idx_dz_dt)

                    % some miscallaneous fields have non-timed information
                    % keep all this information
                    vert_profs(profile_num).(fields{ff}) = vert_profs(profile_num).(fields{ff});

                end


            end




        end




    end


end


%% Sometimes multilayered clouds cause the code above to keep two very similar profiles
% Go through all profiles found. If the indexes of one profile are nearly
% identical to another, get rid of one.

% let's set a buffer of 4 indices. If either the starting or ending indices
% are within 4 of one another, delete the shorter profile
buffer_length = 4;
idx2delete = [];

for n1 = 1:profile_num
    for n2 = (n1+1):profile_num

        if (vert_profs(n2).time(1) - vert_profs(n1).time(1)) <= buffer_length ||...
                (vert_profs(n2).time(end) - vert_profs(n1).time(end)) <= buffer_length

            % If this is true, keep the longer profile
            if length(vert_profs(n2).time)>length(vert_profs(n1).time)
                idx2delete = [idx2delete, n1];

            elseif length(vert_profs(n2).time)<length(vert_profs(n1).time)
                idx2delete = [idx2delete, n2];

            else

                % They're the same size! Just pick one (Is this even
                % possible?)
                idx2delete = [idx2delete, n2];

            end


        end


    end
end


% detele indexes associated with profiles nearly identical to another
% profile
vert_profs(idx2delete) = [];


%% If stop_at_max_LWC is true, remove measurements after max LWC for each profile

if stop_at_max_lwc == true

    error([newline, 'The code to cut profiles after max LWC hasnt been tested', newline])
    

    for nn = 1:length(vert_profs)

        % find index with max LWC
        [max_lwc, idx_max_lwc] = max(vert_profs(nn).lwc);


        for ff = 1:length(fields)

                % now remove all data points outside of the vertical profile

                if numel(vert_profs(nn).(fields{ff}))==length(idx_dz_dt)

                    % all the time data will have be a vector with the same
                    % length as our calculation of vertical velocity

                    % reshape all fields so that time increases with increasing
                    % column number (row vector)
                    vert_profs(nn).(fields{ff}) = vert_profs(nn).(fields{ff})(1:idx_max_lwc);


                elseif numel(vert_profs(nn).(fields{ff}))>length(idx_dz_dt)
                    % only one field is a matrix with more values than the time
                    % vector, and thats the matrix for the size distribution,
                    % with rows representing different size bins and columns
                    % are along the time dimension
                    vert_profs(nn).(fields{ff}) = vert_profs(nn).(fields{ff})(:, 1:idx_max_lwc);


                end


        end


    end

end



%%
for nn = 1:length(vert_profs)

    % Store the LWC threshold used
    vert_profs(nn).LWC_threshold = LWC_threshold;            % g/m^3

    % Store the Nc threshold used
    vert_profs(nn).Nc_threshold = Nc_threshold;            % g/m^3

end



%%
% ----------------------------------------------------------------------
% ------------------ Compute Liquid Water Path -------------------------
% ----------------------------------------------------------------------



% we want to compute the total LWP, and the LWP for the two instruments
% used to measure droplets.
% The CDP probe measures radii between 0.875 and 25.055 microns.
% The 2DC probe measures radii between 31.25 and 793 microns

% step through each profile
for nn = 1:length(vert_profs)


    % LWP is calculated by integrating from cloud bottom to
    % cloud top. If the plane decreases in altitude, we need to
    % integrate from the end of the profile (cloud bottom) to the
    % begining (cloud top). If the plane is ascending, we do
    % the opposite.

    dz_dt = diff(vert_profs(nn).altitude,1)./diff(vert_profs(nn).time, 1);

    if mean(dz_dt)>0
        % then the plane is ascending!

        % Compute the total LWP
        vert_profs(nn).lwp = trapz(vert_profs(nn).altitude, vert_profs(nn).lwc);            % g/m^2

        % ------ Compute the CDP LWP ---------
        vert_profs(nn).lwp_CDP = trapz(vert_profs(nn).altitude, vert_profs(nn).lwc_CDP);


        % ------ Compute the 2DC LWP ---------
        vert_profs(nn).lwp_2DC = trapz(vert_profs(nn).altitude, vert_profs(nn).lwc_2DC);


    elseif mean(dz_dt)<0
        % then the plane is descending!

        % Compute the total LWP
        vert_profs(nn).lwp = trapz(flipud(vert_profs(nn).altitude), flipud(vert_profs(nn).lwc));            % g/m^2

        % ------ Compute the CDP LWP ---------
        vert_profs(nn).lwp_CDP = trapz(flipud(vert_profs(nn).altitude), flipud(vert_profs(nn).lwc_CDP));


        % ------ Compute the 2DC LWP ---------
        % compute the 2DC LWP by integration over the cloud depth
        vert_profs(nn).lwp_2DC = trapz(flipud(vert_profs(nn).altitude), flipud(vert_profs(nn).lwc_2DC));





    end



end








% ----------------------------------------------------------------------
% ----------- COMPUTE THE HORIZONTAL DISTANCE TRAVELLED ----------------
% ----------------------------------------------------------------------

% Also compute the cloud depth and then the slant path travelled in the
% cloud by the plane us pythagreous' theorem. The compute the zenith angle
% of the slant path with repsect to the cloud base.


% why do we need this? Optical depth is a path dependent quantity. The
% properties of the medium that make it attenuating must be integrated
% along the entire line of sight. Cloud optical depth is defind as the
% vertically integrated optical depth from cloud top to cloud bottom.
% Assuming a plane-parallel cloud, the equation is:
% int( sec(theta) * Qe * pi * r^2(z) * Nc(z)) dz
% The term sec(theta) accounts for the zenith angle off from a true
% vertical path (one that travels against the vector of the gravitational
% force). To estimate what sec(theta) is, we need to compute the horizontal
% distance travelled during the vertical profile. Cos(theta) is the ratio
% of the the vertical distance travelled to the hypotenuse. And sec(theta)
% is 1/cos(theta).


% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");        % units of meters


for nn = 1:length(vert_profs)

    % Step through each point to get the linear distance travelled as a
    % vector
    horz_distance_travelled = zeros(1, length(vert_profs(nn).latitude));

    for xx = 2:length(vert_profs(nn).latitude)

        % Find the linear distance between the start and end of the horizontal profile.
        % When you specify a reference ellipsoid as input to the distance function,
        % the function returns linear distances in the units of the semimajor axis of the ellipsoid.
        horz_distance_travelled(xx) = distance(vert_profs(nn).latitude(1), vert_profs(nn).longitude(1),...
            vert_profs(nn).latitude(xx), vert_profs(nn).longitude(xx),wgs84);

    end

    vert_profs(nn).horz_dist = horz_distance_travelled;

    % Using the horizontal distance travelled, and the vertical depth of
    % the cloud, use pythagreous' theorem to estimate the slant path
    % travelled within the cloud
    vert_profs(nn).cloud_depth = max(vert_profs(nn).altitude) - min(vert_profs(nn).altitude);
    vert_profs(nn).slant_path = sqrt(vert_profs(nn).horz_dist(end)^2 + vert_profs(nn).cloud_depth^2);

    % Compute the zenith angle of the slant path with respect to the cloud
    % base, assuming a plane-parallel cloud and a straight line for the
    % planes trajectory

    vert_profs(nn).VR_zenith_angle = atand(vert_profs(nn).horz_dist(end)/vert_profs(nn).cloud_depth);


end




% ------------------------------------------------------------------
% ------------------ Compute optical depth -------------------------
% ------------------------------------------------------------------

% compute the optical depth for each vertical profile and store it


% optical depth is defined to be 0 at cloud top and increasing towards
% cloud bottom

% We do NOT include the zenith angle of the plane as it travels through the
% cloud because we are concerned with the cloud optical depth and not the
% optical depth along the slant path. Cloud optical depth is defined soley
% as the optical depth over the vertical dimension of a cloud.



% step through each profile
for nn = 1:length(vert_profs)


    vector_length = length(vert_profs(nn).altitude);
    vert_profs(nn).tau = zeros(1,vector_length-1);

    % compute sec(theta) by first computing the hypotenuse



    % Optical thickness is defined by integrating from cloud top to
    % cloud bottom. If the plane increases in altitude, we need to
    % integrating from the end of the profile (cloud top) to the
    % begining (cloud bottom). If the plane is descending, we do
    % the opposite.

    dz_dt = diff(vert_profs(nn).altitude)./diff(vert_profs(nn).time);


    if mean(dz_dt)>0
        % then the plane is ascending!


        % step through the altitude array
        for ii = 1:vector_length-1


            % we have to convert Nc and re to have the same units as the alitude,
            % which is in meters

            if vocalsRex.flag_2DC_data_is_conforming==true
                re_meters = vert_profs(nn).re(vector_length-ii:vector_length)./1e6;                      % meters
            else
                % What choice do we have? I guess we will just use the
                % effevtive radius from the 2DC data, but this will
                % underestimate the optical depth
                % if the re CDP data is taken at 10 samples per second,
                % only take every 10th value.
                if length(vocalsRex.re_CDP)>length(vocalsRex.time)
                    error([newline, 'I dont know how to handle SPS10 data!', newline])

                else
                    re_meters = vert_profs(nn).re_CDP(vector_length-ii:vector_length)./1e6;                      % meters
                end

            end

            total_Nc_meters = vert_profs(nn).total_Nc(vector_length-ii:vector_length).*1e6;                           % #/m^3
            altitude = vert_profs(nn).altitude(end) -  vert_profs(nn).altitude(vector_length-ii:vector_length); % meters


            % We assume the droplet size is appreciably larger than the
            % incident wavelength (something in the visible, like 550 nm)
            % so we can assume the extinction efficiency is 2. This leads
            % to the following equation for optical depth:
            % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
            % where Q_e is set to 2
            %vert_profs(nn).tau(ii) = 2*pi* trapz(flipud(altitude), flipud(re_meters.^2 .* total_Nc_meters));

            % Or we can use a pre-computed mie table to use a more
            % accurate value for the extinction efficiency.
            %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters.*1e6], 'mono', true);

            % Or we could compute the average extinction efficiency
            % over a droplet size distrubution
            [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
                550, 'water', 'gamma', which_computer, ii);

            vert_profs(nn).tau(ii) = pi* trapz(fliplr(altitude), fliplr(Qe_avg .* re_meters.^2 .* total_Nc_meters));

        end

    elseif mean(dz_dt)<0
        % then the plane is descending!

        % step through the altitude array
        for ii = 1:vector_length-1


            % we have to convert Nc and re to have the same units as the alitude,
            % which is in meters
            if vocalsRex.flag_2DC_data_is_conforming==true
                re_meters = vert_profs(nn).re(1:ii+1)./1e6;                      % meters
            else
                % What choice do we have? I guess we will just use the
                % effevtive radius from the 2DC data, but this will
                % underestimate the optical depth
                re_meters = vert_profs(nn).re_CDP(1:ii+1)./1e6;                      % meters
            end

            total_Nc_meters = vert_profs(nn).total_Nc(1:ii+1).*1e6;                           % #/m^3
            altitude = vert_profs(nn).altitude(1) -  vert_profs(nn).altitude(1:ii+1);   % meters


            % We assume the droplet size is appreciably larger than the
            % incident wavelength (something in the visible, like 550 nm)
            % so we can assume the extinction efficiency is 2. This leads
            % to the following equation for optical depth:
            % tau = integral( Q_e * pi * r^2(z) * N_c(z) ) dz
            % where Q_e is set to 2
            %vert_profs(nn).tau(ii) = 2*pi* trapz(altitude, re_meters.^2 .* total_Nc_meters);


            % Or we can use a pre-computed mie table to use a more
            % accurate value for the extinction efficiency.
            %Q_e = interp_mie_computed_tables([linspace(550, 550, length(re_meters))', re_meters'.*1e6], 'mono', true);


            % Or we could compute the average extinction efficiency
            % over a droplet size distrubution
            [~, Qe_avg, ~] = average_mie_over_size_distribution(re_meters.*1e6, linspace(10,10,length(re_meters)),...
                550, 'water', 'gamma', which_computer, ii);


            vert_profs(nn).tau(ii) = pi* trapz(altitude, Qe_avg(:,end) .* re_meters.^2 .* total_Nc_meters);



        end

    end



    % add a zero at the begining!
    vert_profs(nn).tau = [0,vert_profs(nn).tau];

end















end