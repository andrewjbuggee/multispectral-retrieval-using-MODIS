function GN_inputs = create_gauss_newton_inputs(modisInputs)

% ----- Define the number of pixels to estimate a profile for -----
%bayes_inputs.numPixels2Calculate = 4;
% For now, let's only compute this calculation for the MODIS pixel in
% question
% We will only use what is in the truth table!
GN_inputs.numPixels2Calculate = modisInputs.pixels.num_2calculate;
% -------------------------------------------------------------------


% define the number of iterations for the gauss-newton solver
GN_inputs.GN_iterations = 5;


% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
GN_inputs.percent_change_limit = 0.03;

% define the type of model prior pdf
GN_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
GN_inputs.num_model_parameters = 3;


% Define the spectral channels to use in the gauss-newton solver

GN_inputs.bands2use = 1:7;  % number of spectral bands to use




% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

% Using the King and Vaughn (2012) method, we retireve 3 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value



    

GN_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth'};


% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


GN_inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
GN_inputs.measurement.covariance_type = 'independent';

% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

% we have to assume a vertical profile of droplet size with cloud optical
% depth exsists. And we only retrieve the droplet size at the top and
% bottom. This is the method of King and Vaughn (2012)

% the options for the vertical droplet profile are:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.
% x is determined by the choice of droplet profile within the function
% create_droplet_profile.m

GN_inputs.model.profile.type = 'adiabatic';
GN_inputs.model.profile.r_top = 10; % microns - value for our model
GN_inputs.model.profile.r_bottom = 5; % microns - value for our model



% -----------------------------------------------
% ------------- Folder Locations  ---------------
% -----------------------------------------------
GN_inputs.save_calcs_fileName = ['uvspec_CALCS_4Bayes_',date,'.mat'];




% -----------------------------------------------
% --------------- Define flags  -----------------
% -----------------------------------------------










% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------


% Define the number of streams to use in your radiative transfer model
GN_inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
GN_inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
GN_inputs.RT.band_parameterization = 'reptran coarse';
%band_parameterization = 'reptran_channel modis_terra_b07';
% ------------------------------------------------------------------------


% ---------------------------------------------------------
% ------ Define the Solar Flux file and it's resolution ---
% ---------------------------------------------------------

% Define the source file
%bayes_inputs.RT.source_file = '../data/solar_flux/kurudz_1.0nm.dat';
% resolution should match the value listed in the file name
% bayes_inputs.RT.sourceFile_resolution = 1;                  % nm

% this is a hybrind reference spectrum downloaded from LASP's
% LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
% These data range from 202 to 2730 nm
% These data have 0.1 sampling resolution
GN_inputs.RT.sourceFile = '../data/solar_flux/hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
% resolution should match the value listed in the file name
GN_inputs.RT.sourceFile_resolution = 0.1;                  % nm

% define the atmospheric data file
GN_inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
GN_inputs.RT.surface_albedo = 0.05;

% day of the year
GN_inputs.RT.day_of_year = str2double(modisInputs.L1B_filename(15:17));




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
GN_inputs.RT.yesCloud = true;

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
GN_inputs.RT.use_MODIS_cloudTopHeight = false;

% --- Do you want to use the VOCALS-REx cloud top height measurement? ----
GN_inputs.RT.use_VOCALS_cloudTopHeight = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS above cloud water vapor? ---------
GN_inputs.RT.use_MODIS_aboveCloudWaterVapor = false;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------
GN_inputs.RT.modify_CO2 = true;

GN_inputs.RT.CO2_mixing_ratio = 416;       % ppm
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------- Do you want to use the VOCALS measured cloud depth? -----------
GN_inputs.RT.use_VOCALS_cloudDepth = true;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
GN_inputs.RT.use_custom_mie_calcs = false;
% ------------------------------------------------------------------------
% This string is used to compute the LWC from optical depth and effective radius
% can be 'hu' or 'mie interpolate'
GN_inputs.RT.wc_parameterization = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% define the type of droplet distribution
GN_inputs.RT.drop_distribution_str = 'gamma';
% define the distribution varaince
% 7 is the value libRadTran uses for liquid water clouds
GN_inputs.RT.drop_distribution_var = 10;
% define whether this is a vertically homogenous cloud or not
GN_inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';
% define how liquid water content will be computed
% can either be 'mie' or '2limit'
GN_inputs.RT.parameterization_str = 'mie';     % This string is used to compute the LWC from optical depth and effective radius


% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
GN_inputs.RT.use_coxMunk = true;
GN_inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
GN_inputs.RT.yesAerosols = true;

GN_inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
GN_inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------


% ----- Do you want a long error message? -----
% if so, set error message to 'verbose'. Otherwise, set error message to
% 'quiet'
GN_inputs.RT.err_msg = 'quiet';










end