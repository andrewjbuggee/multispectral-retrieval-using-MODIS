%% Read in Vocals Rex Data


% By Andrew John Buggee
%%

function vocalsRex = readVocalsRex(filename)

info = ncinfo(filename);

% read the time vector
time = double(ncread(filename, 'Time'));                   % Measured in seconds

% make sure the time data is a row vector
time = reshape(time, 1, length(time));

% [hours, minutes]
startTime = [str2double(info.Variables(1).Attributes(3).Value(26:27)),...
    str2double(info.Variables(1).Attributes(3).Value(29:30))];      % Start time of the data log in UTC


% ----------------------------------------------------------------
% ----------- cloud droplet probe measurements (CDP) -------------
% ----------------------------------------------------------------

% What we measure IS the number concentration for each bin! we do NOT measure the differential dN/dr

num_concentration_CDP = ncread(filename, 'CCDP_RWO');          % #/cm^3 - number concentration for each bin

% The droplet number concentration is filtered into 30 bins from 2 to 52
% microns
% Below 20 microns, the bins are spaced by 1.18 microns
% above 22 microns, the bins are spaced by 2.28 microns
% find the variable information for the droplet concentration data
CCDP_RWO_info = ncinfo(filename, 'CCDP_RWO');

% ----- Check that the data quality is labeled 'Good' -----
if strcmp(CCDP_RWO_info.Attributes(6).Value, 'Good')~= true
    error([newline, 'CDP data quality listed as not good!', newline])
end
% -------

% Bins are defined as diameters. Divide by 2 to get the radius
drop_radius_bin_edges_CDP = double(CCDP_RWO_info.Attributes(10).Value./2);                        % microns

% The data tells us which bins to take
first_bin_CDP = CCDP_RWO_info.Attributes(8).Value;
last_bin_CDP = CCDP_RWO_info.Attributes(9).Value;

% There are 2 edges for each measurement
%
drop_radius_bin_edges_CDP = drop_radius_bin_edges_CDP(first_bin_CDP:last_bin_CDP+1);                         % microns
drop_radius_bin_center_CDP = drop_radius_bin_edges_CDP(1:end-1) + diff(drop_radius_bin_edges_CDP)/2;         % microns

%drop_size_dist_CDP2 = drop_size_dist_CDP;
num_concentration_CDP = num_concentration_CDP(first_bin_CDP:last_bin_CDP,1,:);          % #/cm^3




% ----------------------------------------------------------------
% ---------- Two-Dimensional Optical array probe (2DC) -----------
% ----------------------------------------------------------------

% This droplet probe filters droplets into 64 different size bins. The bins
% range from 25 microns to 1560 microns.

% ** IMPORTANT ** The first bin of the 2DC probe overlaps with the last bin
% of the CDP. So we ignore it. (See Painemal and Zuidema paragraph 15.)

% There have been some issues with naming the 2DC data different variables.
% Let's loop through all the variables and make sure we grab the correct
% one

for nn = 1:length(info.Variables)

    if strcmp('C1DCA_RPC', info.Variables(nn).Name)==true

        num_concentration_2DC = ncread(filename, 'C1DCA_RPC');        % #/L  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!

        % convert the 2DC data to inverse cubic centimeters
        num_concentration_2DC = num_concentration_2DC ./ 1000;         % # of droplets/cm^3 per bin

        % read in the 2DC variable info
        C1DCA_RPC_info = ncinfo(filename, 'C1DCA_RPC');

        % Check to see if the data quality is labeled as 'Good'
        if strcmp(C1DCA_RPC_info.Attributes(6).Value, 'Good')~=true

            error([newline, 'The 2DC data quality is not listed as Good.', newline])

        end


        % Grab the droplet size bins
        % These are the boundaries that define each bin size
        % Bins are defined as diameters. Divide by 2 to get the radius
        drop_radius_bin_edges_2DC = C1DCA_RPC_info.Attributes(12).Value./2;                        % microns

        % The data tells us which bins to take
        for jj = 1:length(C1DCA_RPC_info.Attributes)

            if strcmp(C1DCA_RPC_info.Attributes(jj).Name, 'FirstBin')==true
                first_bin_2DC = double(C1DCA_RPC_info.Attributes(jj).Value);
            end

            if strcmp(C1DCA_RPC_info.Attributes(jj).Name, 'LastBin')==true
                last_bin_2DC = double(C1DCA_RPC_info.Attributes(jj).Value);
            end


        end

        % Sometimes there isn't a frist or last bin. In that case, set them
        % to be the full breadth of the data
        if exist("first_bin_2DC", 'var')==false
            % set it to be 1
            first_bin_2DC = 1;
        end

        if exist("last_bin_2DC", 'var')==false
            % set it to be the length of drop_bin_edges
            last_bin_2DC = length(drop_radius_bin_edges_2DC)-1;
        end


        drop_radius_bin_edges_2DC = drop_radius_bin_edges_2DC(first_bin_2DC:last_bin_2DC+1);                          % microns
        drop_radius_bin_center_2DC = drop_radius_bin_edges_2DC(1:end-1) + diff(drop_radius_bin_edges_2DC)/2;          % microns

        num_concentration_2DC = num_concentration_2DC(first_bin_2DC:last_bin_2DC,1,:);



    end


end


% ------------------------- IMPORTANT ----------------------------------

% if after the entire for loop this variable is not found, we will look for
% another variable containg 2DC data. For some reason, certain data sets
% don't have the variable listed above. Instead that usually have a
% variable described as '2D-C Concentration, 260X Emulation (per cell)'.
% But even this variable is sometimes empty! It's all rather confusing.
% What I have found is that when the above 2DC variable is missing, the
% variable C1DC_RPC often exists and the dependent variables of LWC and
% mean droplet diameter are often populated with values. Let's download
% these variables if the above variable is not found.



if exist("num_concentration_2DC", 'var')==false


    for nn = 1:length(info.Variables)


        if strcmp('C1DC_RPC', info.Variables(nn).Name)==true


            % read in the 2DC variable info
            C1DC_RPC_info = ncinfo(filename, 'C1DC_RPC');
            % First check to make use the data quality is listed as 'Good'.
            % If it is, then proceed

            if strcmp(C1DC_RPC_info.Attributes(6).Value, 'Good')==true

                % If the data quality is good, let's read in the data
                % Read in the C1DC_RPC data

                % set up a flag because we will have to calcualte things
                % differently
                flag_2DC_data_is_conforming = false;


                %raw_counts_2DC = ncread(filename, 'A1DC_RPC');        % counts/bin  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!
                num_concentration_2DC = ncread(filename, 'C1DC_RPC');        % #/L  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!


                % Grab the droplet size bins
                % These are the boundaries that define each bin size
                % Bins are defined as diameters. Divide by 2 to get the radius
                drop_radius_bin_edges_2DC = C1DC_RPC_info.Attributes(12).Value./2;                        % microns

                % The data tells us which bins to take
                first_bin_2DC = double(C1DC_RPC_info.Attributes(10).Value);
                last_bin_2DC = double(C1DC_RPC_info.Attributes(11).Value);

                drop_radius_bin_edges_2DC = drop_radius_bin_edges_2DC(first_bin_2DC:last_bin_2DC+1);                          % microns
                drop_radius_bin_center_2DC = drop_radius_bin_edges_2DC(1:end-1) + diff(drop_radius_bin_edges_2DC)/2;          % microns


                % --- define the 2DC droplet probe data ----
                num_concentration_2DC = num_concentration_2DC(first_bin_2DC:last_bin_2DC,1,:);

                % ---- check to ensure the 2DC data is nan free ----
                if any(isnan(num_concentration_2DC), 'all')==true

                    % then we have no choice but to set the 2DC data to be
                    % a bunch of zeros
                    num_concentration_2DC = zeros(size(num_concentration_2DC));


                end


                % --- Let's also read the LWC data and the mean diameter data
                lwc_2DC = reshape(ncread(filename, 'PLWC1DC_RPC'), 1, []);             % g/m^3 - 2DC liquid water content
                mean_radius_2DC = reshape(ncread(filename, 'DBAR1DC_RPC'), 1, [])./2;   % microns - 2DC mean particle radius

                % we should alsos read in the total droplet concentration data
                total_Nc_2DC = ncread(filename, 'CONC1DC_RPC');     % #/L
                total_Nc_2DC = reshape(total_Nc_2DC./ 1000, 1, []);                 % # of droplets/cm^3


            end


        end


    end

else

    % if this is true, then we have good 2DC data to work with. Tell the code
    % the 2DC data is conforming.

    % set up a flag because we will have to calcualte things
    % differently
    flag_2DC_data_is_conforming = true;


end



% Let's check one more time. There are two 2DC probes that we could use. If
% the above probe did not exist in the file, or if the data qaulity was
% bad, we can try one more...

if exist("num_concentration_2DC", 'var')==false


    for nn = 1:length(info.Variables)


        if strcmp('C1DC_RPI', info.Variables(nn).Name)==true


            % read in the 2DC variable info
            C1DC_RPI_info = ncinfo(filename, 'C1DC_RPI');
            % First check to make use the data quality is listed as 'Good'.
            % If it is, then proceed

            if strcmp(C1DC_RPI_info.Attributes(6).Value, 'Good')==true

                % If the data quality is good, let's read in the data
                % Read in the C1DC_RPC data

                % set up a flag because we will have to calcualte things
                % differently
                flag_2DC_data_is_conforming = false;


                %raw_counts_2DC = ncread(filename, 'A1DC_RPI');        % counts/bin  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!
                num_concentration_2DC = ncread(filename, 'C1DC_RPI');        % #/L  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!


                % Grab the droplet size bins
                % These are the boundaries that define each bin size
                % Bins are defined as diameters. Divide by 2 to get the radius
                drop_radius_bin_edges_2DC = C1DC_RPI_info.Attributes(12).Value./2;                        % microns

                % The data tells us which bins to take
                first_bin_2DC = double(C1DC_RPI_info.Attributes(10).Value);
                last_bin_2DC = double(C1DC_RPI_info.Attributes(11).Value);

                drop_radius_bin_edges_2DC = drop_radius_bin_edges_2DC(first_bin_2DC:last_bin_2DC+1);                          % microns
                drop_radius_bin_center_2DC = drop_radius_bin_edges_2DC(1:end-1) + diff(drop_radius_bin_edges_2DC)/2;          % microns


                % --- define the 2DC droplet probe data ----
                num_concentration_2DC = num_concentration_2DC(first_bin_2DC:last_bin_2DC,1,:);


                % ---- check to ensure the 2DC data is nan free ----
                if any(isnan(num_concentration_2DC), 'all')==true

                    % then we have no choice but to set the 2DC data to be
                    % a bunch of zeros
                    num_concentration_2DC = zeros(size(num_concentration_2DC));


                end


                % ----------------- Read the LWC data --------------
                lwc_2DC = reshape(ncread(filename, 'PLWC1DC_RPI'), 1, []);              % g/m^3 - 2DC liquid water content

                % -------- We cannot compute the effective radius --------
                % Wihtout number concentration data, we cannot compute the
                % 2DC effective radius. And I haven't found it in the list
                % of variables. But the do compute a 'mean' radius, which
                % appears to be very close to the first moment of the
                % droplet distribution: int(r*n(r))dr / int(n(r))dr
                mean_radius_2DC = reshape(ncread(filename, 'DBAR1DC_RPI'), 1, [])./2;      % microns - 2DC mean particle radius

                % we should alsos read in the total droplet concentration data
                total_Nc_2DC = ncread(filename, 'CONC1DC_RPI');     % #/L
                total_Nc_2DC = reshape(total_Nc_2DC./ 1000, 1, []);           % # of droplets/cm^3

            end


        end


    end


end


% If we STILL don't have any 2DC data, throw an error

if exist("num_concentration_2DC", 'var')==false

    error([newline, 'There is not 2DC data available in this dataset.', newline])


end




% -------------------------------------------------------------------------
% ------- Grab Wind speed, direction and sensor lat,long, altititude ------
% -------------------------------------------------------------------------

% grab the horizontal wind speed and direction
horz_wind_speed = ncread(filename, 'WSC');                                    % m/s -  GPS corrected horizontal wind speed
horz_wind_direction = ncread(filename, 'WDC');                                 % degrees from north -  GPS corrected horizontal wind direction 'wind from'

% Grab position and altitude data
lat = ncread(filename, 'LAT');                                                        % Meausred in degrees North
long = ncread(filename, 'LON');                                                       % Meausred in degrees East
altitude = ncread(filename, 'ALTX');                                                  % Measured in meters above sea level


% Somestimes the altitude vector is not a vector at all but a matrix.
% Usually each row is roughly the same. Just take the first row
if numel(altitude)~=numel(time)
    altitude = altitude(1,:);
end

% Make sure that altitude is a row vector
altitude = reshape(altitude, 1, length(altitude));


% -------------------------------------------------------------------------
% ------------------ Grab the ambient air temperature ---------------------
% -------------------------------------------------------------------------
ambient_air_temp = reshape(ncread(filename, 'ATX'), 1, []);                                          % C - ambient air tempertaure


% -------------------------------------------------------------------------
% --------------- Grab different water vapor measurements ---------------
% -------------------------------------------------------------------------
water_vapor_pressure = reshape(ncread(filename, 'EDPUV'), 1, []);                                % hPa - ambient water vapor pressure
mixing_ratio = reshape(ncread(filename, 'MR'), 1, []);                                           % g of water vapor/kg of dry air - water vapor mixing ratio
absolute_humidity = reshape(ncread(filename, 'RHODT'), 1, []);                                   % g of water vapor/m^3 of air - absolute humidity





% -------------------------------------------------------------------------
% ------- Grab the short and longwave radiances looking up and down ------
% -------------------------------------------------------------------------

% lets grab the shortwave and longwave radiance looking down and looking up
shortwave_top = ncread(filename, 'SWT');
shortwave_bot = ncread(filename, 'SWB');

longwave_top = ncread(filename,'IRTC');              % These are corrected longwave irradiance values
longwave_bot = ncread(filename, 'IRBC');






% -------------------------------------------------------------------
% ----------------- ***** IMPORTANT STEP ***** ----------------------
% -------------------------------------------------------------------


% Combine data from both droplet probes to get the total droplet size
% distribution

% ***** There is a bin where no data exists. The CDP data ends at bin
% [24.0550 - 25.195] microns. The 2DC data starts with the bin
% [31.25 - 43.75] microns. Therefore we have to place a 0 inbetween these
% two data bins. That is, since we have no measurement, we define the
% number concentration to be 0 between [25.195 - 31.25] *****

% we also need to define the center point of the bin that has 0
% measurements
center_for_0 = drop_radius_bin_edges_CDP(end) + (drop_radius_bin_edges_2DC(1) - drop_radius_bin_edges_CDP(end))/2;                               % microns
drop_radius_bin_center = [drop_radius_bin_center_CDP, center_for_0, drop_radius_bin_center_2DC];                    % microns                                                   % microns
drop_radius_bin_edges = [drop_radius_bin_edges_CDP, drop_radius_bin_edges_2DC];                                                             % microns
%drop_radius_bin_edges2 = [drop_radius_bin_edges_CDP2, drop_radius_bin_edges_2DC2];                                                             % microns
%drop_radius_bin_center2 = drop_radius_bin_edges2(1:end-1) + diff(drop_radius_bin_edges2)./2;                                                  % microns
%drop_radius_bin_first_edge = [drop_radius_bin_edges_CDP(1:end-1), drop_radius_bin_edges_2DC(1:end-1)];                                      % microns

% don't forget the 0!
% Nc is measured in #/cm^3
Nc = [reshape(num_concentration_CDP,length(drop_radius_bin_center_CDP),[]); zeros(1,length(time)); reshape(num_concentration_2DC,length(drop_radius_bin_center_2DC),[])];          % Number of droplets in each bin
%Nc2 = [reshape(drop_size_dist_CDP2,length(drop_radius_bin_edges_CDP2),[]); reshape(drop_size_dist_2DC2,length(drop_radius_bin_edges_2DC2),[])];          % Number of droplets in each bin

droplet_matrix_center = repmat((drop_radius_bin_center)', 1, length(time))./1e4;            % cm                                                         % cm
%droplet_matrix_center2 = repmat((drop_radius_bin_center2)', 1, length(time));                                                                  % microns
%droplet_matrix_firstEdge = repmat((drop_radius_bin_first_edge)', 1, length(time))./1e4;                                                        % cm
%droplet_matrix_edges = repmat((drop_radius_bin_edges)', 1, length(time))./1e4;                                                                       % microns
%droplet_matrix_edges2 = repmat((drop_radius_bin_edges2)', 1, length(time))./1e4;                                                                       % microns


% -------------------------------------------------------------------
% ----------- To compute efffective radius we need dN/dr ------------
% -------------------------------------------------------------------

%dN_dr = diff(Nc2,1,1)./diff(droplet_matrix_edges2,1,1);                       % #/micron/cm^3     this is the droplet distribution

% Lets compute the effective radius, which is the 3rd moment to the 2nd
% moment

if flag_2DC_data_is_conforming==true

    % --- IF TRUE, WE HAVE REAL 2DC DATA, NOT THE 260X EMULATION DATA ---



    % ------------------------------------------------------------------
    % ---------- Compute the total Number Concentration ---------------
    % ------------------------------------------------------------------
    % Lets compute the total number concetration at each time step by
    % integrating over r
    %total_Nc = trapz(drop_bins, Nc,1);       % cm^(-3)
    total_Nc = double(sum(Nc,1));                     % cm^(-3)

    % Now compute the effective radius for just the CDP instrument
    index_r_cdp = (drop_radius_bin_center<=drop_radius_bin_center_CDP(end))';       % preform this in microns

    % Now compute the total number concentration for the CDP instrument
    total_Nc_CDP = double(sum(Nc(index_r_cdp,:),1));                     % cm^(-3)

    % Now compute the total number concentration for the 2DC instrument
    total_Nc_2DC = double(sum(Nc(~index_r_cdp,:),1));                     % cm^(-3)



    % ------------------------------------------------------------------
    % --------------- Compute liquid water content ---------------------
    % ------------------------------------------------------------------

    % Lets compute the total liquid water content at each time

    % Lets compute the liquid water content and liquid water path
    rho_lw = 1e6;                                                   % g/m^3 - density of liquid water

    % we have to convert re to cm in order to have the finals units be in grams
    % per meter cubed

    lwc = double( 4/3 * pi *  rho_lw * sum(Nc .* droplet_matrix_center.^3,1) );        % grams of liquid water/meter cubed of air

    % **** Painemal and Zuidema (2011) found a bias in the CDP flight
    % measurements of LWC. To correct for this, they set the CDP LWC values
    % to be equal to the measurements made by the King hot wire probe,
    % which they claim has a high correlation with the Gerber PV-100 probe.
    % We will preform a similar correction by setting the CDP LWC to be
    % equal to the PV-100 probe data, which is consistantly less than the
    % CDP measurements.
    lwc_PV100 = reshape(ncread(filename, 'XGLWC'), 1, []);                         % g/m^3 - Gerber PV-100 Probe Liquid Water Content

    % Set values less than 0 to be 0
    lwc_PV100(lwc_PV100<0) = 0;
    % ------------------------- CDP LWC ----------------------------
    % compute the liquid water content measured by the CDP Instrument
    lwc_CDP = double( 4/3 * pi *  rho_lw * sum(Nc(index_r_cdp, :) .* droplet_matrix_center(index_r_cdp, :).^3,1) );     % grams of liquid water/meter cubed of air

    % Solve for the coefficient (Painemal and Zuidema 2011 pg 4)
    a = lwc_CDP./lwc_PV100;

    % set NaN values to 1
    a(isnan(a)) = 1;

    % set the inf calues to 1
    a(a==inf) = 1;

    % set zero values to be 1
    a(a==0) = 1;

    % Compute the corrected LWC values for the CDP instrument
    % According to Painemal and Zuidema 2011 pg 4, use a to correct the LWC
    % bias by creating a modified center radius r' = (r/a^(1/3))
    lwc_CDP = double( 4/3 * pi *  rho_lw * sum(Nc(index_r_cdp, :) .* (droplet_matrix_center(index_r_cdp, :)./a.^(1/3)).^3,1) );     % grams of liquid water/meter cubed of air




    % ------------------------- 2DC LWC ----------------------------
    % compute the liquid water content measured by the 2DC Instrument
    lwc_2DC = double( 4/3 * pi *  rho_lw * sum(Nc(~index_r_cdp, :) .* droplet_matrix_center(~index_r_cdp, :).^3,1) );                  % grams of liquid water/meter cubed of air

    %     if exist('C1DC_RPC_info', 'var')==true || exist('C1DCA_RPC_info', 'var')==true
    %         lwc_2DC = reshape(ncread(filename, 'PLWC1DC_RPC'), 1, []);
    %
    %     elseif exist('C1DC_RPI_info', 'var')==true
    %         lwc_2DC = reshape(ncread(filename, 'PLWC1DC_RPI'), 1, []);
    %     end



    % *** Now recalculate the total LWC **

    lwc = lwc_CDP + lwc_2DC;            % g/m^3





    % ------------------------------------------------------------------
    % ----------------- Compute the Effective Radius -------------------
    % ------------------------------------------------------------------

    % Compute the ratio of the third moment to the second moment and convert
    % back to microns

    % Use the correction described by Painemal and Zuidema (2011) page 4.
    % The bias correction comes from solving for the difference between the
    % CDP meausred LWC, which tends to be over estimated, and the King hot
    % wire probe LWC, which in our case is the Gerber PV-100 LWC. The
    % radius correction only applies for the CDP instrument bins.

    droplet_matrix_center(index_r_cdp, :) = droplet_matrix_center(index_r_cdp,:)./...
                                repmat(a.^(1/3), sum(index_r_cdp), 1);                      % cm

    % *** 0 divided by 0 gives NaN. Set these values to zero ****
    re = double(sum(droplet_matrix_center.^3 .* Nc, 1)./sum(droplet_matrix_center.^2 .* Nc,1) * 1e4);                 % microns

    % set NaN values to 0
    re(isnan(re)) = 0;


    % ------------------ Re CDP ---------------------
    % compute the effective radius using only CDP data
    re_CDP = double(sum(droplet_matrix_center(index_r_cdp,:).^3 .* Nc(index_r_cdp, :), 1)./...
        sum(droplet_matrix_center(index_r_cdp, :).^2 .* Nc(index_r_cdp, :),1) * 1e4);                 % microns

    % *** 0 divided by 0 gives NaN. Set these values to zero ****
    % set NaN values to 0
    re_CDP(isnan(re_CDP)) = 0;

    %     if strcmp(filename(end-39:end-35), 'SPS_1')==true
    %             % If we wish to read in 1Hz data, take the median at each time
    %             % step.
    %             re_CDP = ncread(filename, 'REFFD_RWO');         % microns
    %
    %             % Check to make sure we only have 1Hz data. Sometimes we dont!
    %             if size(re_CDP,1)*size(re_CDP,2) == size(time,1)*size(time,2)
    %                 re_CDP = reshape(re_CDP, 1, []);
    %
    %             elseif size(re_CDP,1)*size(re_CDP,2) > size(time,1)*size(time,2)
    %
    %                 % find which dimension has the 10 Hz data
    %                 if size(re_CDP,1)==10
    %                     re_CDP = median(re_CDP, 1);
    %
    %                 elseif size(re_CDP,2)==10
    %                     re_CDP = median(re_CDP, 2);
    %                     re_CDP = reshape(re_CDP, 1, []);
    %                 end
    %
    %             else
    %
    %                 error([newline, 'I dont know what to do with the re_CDP data.', newline])
    %
    %             end
    %
    %
    %         elseif strcmp(filename(40:end-35), 'SPS_25')==true
    %             % If we wish to have 10 Hz data (of which the files are labeled
    %             % SPS 25, then we simply read in all data
    %             re_CDP = reshape(ncread(filename, 'REFFD_RWO'), 1, []);
    %     end

    % ------------------ Re 2DC ---------------------
    % compute the effective radius using only 2DC data
    re_2DC = double(sum(droplet_matrix_center(~index_r_cdp, :).^3 .* Nc(~index_r_cdp, :), 1)./...
        sum(droplet_matrix_center(~index_r_cdp, :).^2 .* Nc(~index_r_cdp, :),1) * 1e4);                 % microns

    % *** 0 divided by 0 gives NaN. Set these values to zero ****
    % set NaN values to 0
    re_2DC(isnan(re_2DC)) = 0;






else

    % If the 2DC is nonconforming, let's first check to see if it has any
    % non-zero values.

    if sum(num_concentration_2DC, "all")==0

        % then the 2DC matrix was defined as an array of zeros.
        % We can compute the effective radius for the 2DC instrument and
        % the CDP instrument seperately, but not together

        % HOWEVER, if the entire 2DC number concentration data set consists
        % of zeros, the computed effective radius value IS the effective
        % radius of the CDP instrument



        % Now compute the effective radius for just the CDP instrument
        index_r_cdp = (drop_radius_bin_center<=drop_radius_bin_center_CDP(end))';       % preform this in microns



        % ------------------------------------------------------------------
        % --------------- Compute liquid water content ---------------------
        % ------------------------------------------------------------------

        % Lets compute the total liquid water content at each time step by
        % first seperating it by instruments
        % Again, because the number  data consists of only
        % zeros for the 2DC instrument, the sum is only the contribution
        % from the CDP instrument

        % **** Painemal and Zuidema (2011) found a bias in the CDP flight
        % measurements of LWC. To correct for this, they set the CDP LWC values
        % to be equal to the measurements made by the King hot wire probe,
        % which they claim has a high correlation with the Gerber PV-100 probe.
        % We will preform a similar correction by setting the CDP LWC to be
        % equal to the PV-100 probe data, which is consistantly less than the
        % CDP measurements.
        lwc_PV100 = reshape(ncread(filename, 'XGLWC'), 1, []);                         % g/m^3 - Gerber PV-100 Probe Liquid Water Content

        % Set values less than 0 to be 0. No negative values allowed
        lwc_PV100(lwc_PV100<0) = 0;
        % ------------------------- CDP LWC ----------------------------

        % Lets compute the liquid water content and liquid water path
        rho_lw = 1e6;                                                   % g/m^3 - density of liquid water

        % compute the liquid water content measured by the CDP Instrument
        lwc_CDP = double( 4/3 * pi *  rho_lw * sum(Nc(index_r_cdp, :) .* droplet_matrix_center(index_r_cdp, :).^3,1) );                  % grams of liquid water/meter cubed of air

        % Solve for the coefficient (Painemal and Zuidema 2011 pg 4
        a = lwc_CDP./lwc_PV100;

        % set NaN values to 1
        a(isnan(a)) = 1;

        % set the inf calues to 1
        a(a==inf) = 1;

        % set zero values to be 1
        a(a==0) = 1;

        % Compute the corrected LWC values for the CDP instrument
        % According to Painemal and Zuidema 2011 pg 4, use a to correct the LWC
        % bias by creating a modified center radius r' = (r/a^(1/3))
        lwc_CDP = double( 4/3 * pi *  rho_lw * sum(Nc(index_r_cdp, :) .* (droplet_matrix_center(index_r_cdp, :)./a.^(1/3)).^3,1) );     % grams of liquid water/meter cubed of air



        % compute the total liquid water content
        % Check to see if the data is SPS1 or SPS10
        if length(lwc_CDP)>length(time)

            lwc = max(ncread(filename, 'PLWCD_RWO'), [], 1) + lwc_2DC;      % g/m^3

        elseif length(lwc_CDP)==length(time)

            lwc = lwc_CDP + lwc_2DC;              % g/m^3

        end


        % Use the correction described by Painemal and Zuidema (2011) page 4.
        % The bias correction comes from solving for the difference between the
        % CDP meausred LWC, which tends to be over estimated, and the King hot
        % wire probe LWC, which in our case is the Gerber PV-100 LWC. The
        % radius correction only applies for the CDP instrument bins.

        droplet_matrix_center(index_r_cdp, :) = droplet_matrix_center(index_r_cdp,:)./...
            repmat(a.^(1/3), sum(index_r_cdp), 1);                      % cm

        % *** 0 divided by 0 gives NaN. Set these values to zero ****
        re = double(sum(droplet_matrix_center.^3 .* Nc, 1)./sum(droplet_matrix_center.^2 .* Nc,1) * 1e4);                 % microns

        % set NaN values to 0
        re(isnan(re)) = 0;


        % ------------------ Re CDP ---------------------
        % compute the effective radius using only CDP data
        re_CDP = double(sum(droplet_matrix_center(index_r_cdp,:).^3 .* Nc(index_r_cdp, :), 1)./...
            sum(droplet_matrix_center(index_r_cdp, :).^2 .* Nc(index_r_cdp, :),1) * 1e4);                 % microns

        % *** 0 divided by 0 gives NaN. Set these values to zero ****
        % set NaN values to 0
        re_CDP(isnan(re_CDP)) = 0;



        % ------------------------------------------------------------------
        % ---------- Compute the total Number Concentration ---------------
        % ------------------------------------------------------------------

        % Lets compute the total number concetration at each time step by
        % first seperating it by instruments
        % Again, because the number concentration data consists of only
        % zeros for the 2DC instrument, the sum is only the contribution
        % from the CDP instrument
        total_Nc_CDP = double(sum(Nc,1));                     % cm^(-3)

        % then the total Nc is the sum of the two
        total_Nc = total_Nc_CDP + total_Nc_2DC;             % # / cm^{-3}









    else

        % Then there is data? And we do it the onld fashion way but we
        % still hold onto the 2DC calculated products.
        error([newline, 'The 2DC data is the weird emulsion kind but there are more than just zeros!', newline])


    end

end



%% Compute water vapor concentration

% open physical constants library
con = physical_constants;

% compute number concentration from absolute humidity
vocalsRex.Nc_vapor = absolute_humidity .* (con.N_A/(con.Mol_mass_h20_vap*1e3));             % # of molecules/m^3

% convert units to inverse cubic cm
vocalsRex.Nc_vapor = vocalsRex.Nc_vapor ./ 1e6;                                             % # of molecules/cm^3

% try computing number concentration from the partial pressure as well
vocalsRex.Nc_vapor_partPres = (water_vapor_pressure.*1e2) .* (con.N_A./(con.R_uni * (ambient_air_temp + 273.15)));     % # of molecules / m^3

% convert units to inverse cubic cm
vocalsRex.Nc_vapor_partPres = vocalsRex.Nc_vapor_partPres ./ 1e6;                            % # of molecules/cm^3




%% Collect all of the data
% We ignore the last bin in the data set. According to the nc info file,
% the data does not extend through all 31 rows.

% store the 2DC flag
vocalsRex.flag_2DC_data_is_conforming = flag_2DC_data_is_conforming;

vocalsRex.Nc = Nc;                                              % both instruments
vocalsRex.drop_radius_bin_edges = drop_radius_bin_edges;
vocalsRex.drop_radius_bin_center = drop_radius_bin_center;
vocalsRex.total_Nc = total_Nc;                                  % both instruments
vocalsRex.lwc = double(lwc);                                    % g/m^3 both instruments

% some 

% we only have droplet effective radius from both instruments if the 2DC
% data is non-zero
if flag_2DC_data_is_conforming==true

    % then we have an effective radius that uses data from both instruments
    vocalsRex.re = re;                                          % both instruments
    % and we have an effective radius from just the 2DC data
    vocalsRex.re_2DC = re_2DC;                                               % from the 2DC instrument only

else

    % we don't have an effective radius for the 2DC data
    % What we we have is the first moment
    vocalsRex.mean_r_2DC = mean_radius_2DC;                                               % from the 2DC instrument only

end
vocalsRex.re_CDP = re_CDP;                                               % from the CDP instrument only
vocalsRex.lwc_CDP = lwc_CDP;                                              % from the CDP instrument only
vocalsRex.total_Nc_CDP = total_Nc_CDP;                                               % from the CDP instrument only
vocalsRex.lwc_2DC = lwc_2DC;                                              % from the 2DC instrument only
vocalsRex.total_Nc_2DC = total_Nc_2DC;                                               % from the 2DC instrument only

% create a time vector for SPS10, if there is need
% if length(lwc_CDP)>length(time)
%     vocalsRex.time = linspace(time(1), time(end), length(lwc_CDP));     % sec
% end
vocalsRex.time = time;                                                  % sec
vocalsRex.startTime = startTime;                                  % We have to assume that this is in UTC time as well
vocalsRex.time_utc = double((startTime(1) + startTime(2)/60) + vocalsRex.time./3600); % UTC time
vocalsRex.latitude = lat;
vocalsRex.longitude = long;
vocalsRex.altitude = altitude;

vocalsRex.horz_wind_speed = horz_wind_speed;
vocalsRex.horz_wind_direction = horz_wind_direction;

vocalsRex.vapor_pressure = water_vapor_pressure;             % hPa - ambient water vapor pressure
vocalsRex.mixing_ratio = mixing_ratio;                       % g of water vapor/kg of dry air - water vapor mixing ratio
vocalsRex.absolute_humidity = absolute_humidity;             % g of water vapor/m^3 of air - absolute humidity

vocalsRex.temp = ambient_air_temp;                           % degrees C


% ------------ STUFF I'M CURERENTLY NOT USING -------------
% All this stuff will be useful someday! But right now, I'm not using it.


% vocalsRex.ambient_air_temp = ambient_air_temp;
% vocalsRex.water_vapor_presure = water_vapor_pressure;
% vocalsRex.SWT = shortwave_top;
% vocalsRex.SWB = shortwave_bot;
% vocalsRex.LWT = longwave_top;
% vocalsRex.LWB = longwave_bot;







end