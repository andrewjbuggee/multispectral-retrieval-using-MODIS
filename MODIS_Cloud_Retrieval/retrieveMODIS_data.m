%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function [modis,L1B_fileNames] = retrieveMODIS_data(folderName)


% --- Add foldername to the matlab path ---
addpath(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName,'*.hdf']); % find all files that end in .hdf

% check to see if we found any files!
if isempty(files)==true
    error([newline,'There are no files in the folder provided!', newline])
end

L1B_fileNames = cell(1,length(files));

for ii = 1:length(files)

    file_ii = files(ii).name;

    if strcmp(file_ii(1:5),'MOD02') == true || strcmp(file_ii(1:5), 'MYD02') == true

        % Record the time the data was taken
        % Time recorded in [hours, minutes]
        modis.time = [str2double(file_ii(19:20)), str2double(file_ii(21:22))];        % UTC time of data recording
        modis.time_decimal = modis.time(1) + modis.time(2)/60;                        % UTC time in decimal format

        if strcmp(file_ii(1:8),'MOD02QKM') == true || strcmp(file_ii(1:8),'MYD02QKM') == true
            % 250m resolution calibrated data

            L1B_fileNames{ii} = file_ii;
            %modis.EV250m = compute_MODIS_sub_pixel_heterogeneity(file_ii);


        elseif strcmp(file_ii(1:8),'MOD02HKM') == true || strcmp(file_ii(1:8),'MYD02HKM') == true
            % 500m resolution calibrated data

            L1B_fileNames{ii} = file_ii;
            modis.EV500m = readMODIS_L1B_data(file_ii);

        elseif strcmp(file_ii(1:8),'MOD021KM') == true || strcmp(file_ii(1:8),'MYD021KM') == true
            % 1km resolution data, which is the same resolution as the
            % cloud microphysical retrievals

            L1B_fileNames{ii} = file_ii;
            modis.EV1km = readMODIS_L1B_data(file_ii);

        end



    elseif strcmp(file_ii(1:5),'MOD03') == true || strcmp(file_ii(1:5), 'MYD03') == true

        % save the filename
        L1B_fileNames{ii} = file_ii;

        % extract geolocation data from hdf files
        [modis.sensor,modis.solar,modis.geo] = readMODIS_geolocation(file_ii);

    elseif strcmp(file_ii(1:5),'MOD06') == true || strcmp(file_ii(1:5), 'MYD06') == true

        % save the filename
        L1B_fileNames{ii} = file_ii;

        % Retrive the true MODIS cloud Properties

        modis.cloud = readMODIS_L2_data(file_ii);


    elseif strcmp(file_ii(1:5), 'MOD05') == true || strcmp(file_ii(1:5), 'MYD05') == true

        % save the filename
        L1B_fileNames{ii} = file_ii;

        % Retrieve the water vapor column estimates
        modis.vapor = readMODIS_L2_waterVapor_data(file_ii);   % values in cm of water precipitable water

    elseif strcmp(file_ii(1:5),'MOD35') == true || strcmp(file_ii(1:5), 'MYD35') == true

        % This is the cloud mask file that contains the computation of H,
        % the horizontal homogeneity index

        % save the filename
        L1B_fileNames{ii} = file_ii;

        modis.cloud.inhomogeneity_index = read_inhomogeneity_index(file_ii);    % percent



    end



end



end