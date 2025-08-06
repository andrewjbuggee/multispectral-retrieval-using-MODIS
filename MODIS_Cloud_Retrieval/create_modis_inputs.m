%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON MODIS DATA -----


% INPUTS:
%   (1) folderName - 

%   (2) L1B_fileName - 


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_modis_inputs(folderName, L1B_fileNames)


% Save the computer name
inputs.which_computer = whatComputer();

% --- SAVE THE MODIS FILE NAME ----
inputs.modisDataFolder = folderName;



% ----- Save the L1B file name -----
inputs.L1B_filename = L1B_fileNames{1};


% ----- defining what pixels to use -----

% we will randomly select this many pixels from the set of suitable pixels
% found to create .INP files Each pixel and its associated geometry will
% have an INP file of its own that will solve the equation of radiative
% transfer
inputs.pixels.num_2calculate = 10;



% only find pixels in the modis data that is greater than or equal to a tau 
% defined by the value below
inputs.pixels.tau_min_threshold = 3; 

% only find pixels in the modis data that is less than or equal to a tau 
% defined by the value below
inputs.pixels.tau_max_threshold = 80; 


% only find pixels in the modis data that have an effective radius greater
% than of equal to the value below
inputs.pixels.re_min_threshold = 3; 


% only find pixels in the modis data that have an effective radius less
% than of equal to the value below
inputs.pixels.re_max_threshold = 24;
 



inputs.bands2run = [1,6,7]; % these are the bands that we will run uvspec with
inputs.bands2search = [1,7; 1,6]; % these are the modis bands that are used in the retrieval problem
inputs.bands2plot = [1,7]; % these are the modis bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.

% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------
inputs.savedCalculations_folderName = [folderName, 'Retrieval_outputs_', char(datetime("today")),'/']; % this is the folder that all the saved calculations will go


% rev_num = 1;
% inputs.savedCalculations_folderName = [folderName, 'Retrieval_outputs_', datetime("today"),...
%     '_rev', num2str(rev_num), '/']; % this is the folder that all the saved calculations will go
% 
% % If this foldername already exists, AND there are .mat files in it, create
% % a new folder by increasing the rev number
% while isfolder(inputs.savedCalculations_folderName)==true
%     
%     rev_num = rev_num +1;
%     inputs.savedCalculations_folderName = [folderName, 'Retrieval_outputs_', datetime("today"),...
%     '_rev', num2str(rev_num), '/']; % this is the folder that all the saved calculations will go
% 
% 
% end

inputs.saveCalculations_fileName = ['uvspec_calculations_', char(datetime("today")),'.mat'];

% save the day of the year
if L1B_fileNames{1}(15)==0
    inputs.day_of_year = L1B_fileNames{1}(16:17);           % day of year of the MODIS measurement

else
    inputs.day_of_year = L1B_fileNames{1}(15:17);           % day of year of the MODIS measurement
end

% Define the folder to save all the INP files in using the month, day and
% year
data_date = datetime([L1B_fileNames{1}(11:14),'-01-01'],'InputFormat','yyyy-MM-dd') + days(str2double(L1B_fileNames{1}(15:17)) -1);
% check to see if the MODIS instrument is aboard Terra or Aqua
if strcmp(L1B_fileNames{1}(1:3), 'MOD')==true
    % Then read in the spectral response functions for the terra instrument
    inputs.INP_folderName = ['MODIS_Terra_',char(data_date),'_time_',L1B_fileNames{1}(19:22),'/']; % this is the folder name that the INP files will be written to 
elseif strcmp(L1B_fileNames{1}(1:3), 'MYD')==true
    % Then read in the spectral response functions for the Aqua instrument
        inputs.INP_folderName = ['MODIS_Aqua_',char(data_date),'_time_',L1B_fileNames{1}(19:22),'/']; % this is the folder name that the INP files will be written to 

end








% ------------------
% ----- FLAGS! -----
% ------------------

% define flags that tell the codes to either run certain things, or don't
% run certain things

inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use

% if true, the code will load an older set of pixels that has already been used before, and 
% likely has INP files. If false, it tells the code to find a new random subset of pixels
inputs.flags.loadPixelSet = false; 
inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot the 





% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------


% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
inputs.RT.band_parameterization = 'reptran coarse';
%band_parameterization = 'reptran_channel modis_terra_b07';
% ------------------------------------------------------------------------


% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------

% Define the source file
%inputs.RT.source_file = '../data/solar_flux/kurudz_1.0nm.dat';
% resolution should match the value listed in the file name
%inputs.RT.sourceFile_resolution = 1;                  % nm

% hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc has 0.1nm
% sampling resolution
inputs.RT.source.file = '../data/solar_flux/hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
% resolution should match the value listed in the file name
inputs.RT.sourceFile_resolution = 0.1;                  % nm


% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.surface_albedo = 0.05;

% day of the year
inputs.RT.day_of_year = str2double(L1B_fileNames{1}(15:17));




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
inputs.RT.linear_cloudFraction = false;
% if false, define the cloud cover percentage
inputs.RT.cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
inputs.RT.use_MODIS_cloudTopHeight = false;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS above cloud water vapor? ---------
inputs.RT.use_MODIS_aboveCloudWaterVapor = false;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
inputs.RT.drop_distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
inputs.RT.drop_distribution_var = 10;
% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius


% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
inputs.RT.yesAerosols = true;

inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
inputs.RT.err_msg = 'quiet';









% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
