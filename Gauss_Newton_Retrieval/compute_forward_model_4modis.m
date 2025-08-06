% ---- Compute fowrad model using MODIS spectral channels ----

% this function will compute reflectances using LibRadTran radiative
% transfer solver. It will do this for a water cloud with a droplet
% profile over MODIS spectral channels

% for now we will run this using the first 7 MODIS channels

% By Andrew J. Buggee
%%
function measurement_estimate = compute_forward_model_4modis(modis,current_guess,GN_inputs,pixel_row,pixel_col,modisInputs, pp)

% Define some needed folder and file names
saveCalculations_fileName = GN_inputs.save_calcs_fileName; % where to save the computed data
INP_folderName = modisInputs.INP_folderName; % Where to save the INP files

% --- compute the forward model at our current estimate ---
r_top = current_guess(1);
r_bottom = current_guess(2);
tau_c = current_guess(3);

profile_type = GN_inputs.model.profile.type; % type of water droplet profile

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = modisBands(1);    % nm - Wavelength used for cloud optical depth calculation
% ----------------------------------------------------------

% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------

% Set up a few constants for the water cloud
H = GN_inputs.RT.cloudDepth;                                % km - geometric thickness of cloud
n_layers = GN_inputs.RT.cloud_layers;                          % number of layers to model within cloud

% Cloud top
z_top = GN_inputs.RT.cloudTop_height(pp);        % km -  cloud top height

z = linspace(z_top-H, z_top,n_layers);        % km - altitude above ground vector

indVar = 'altitude';                    % string that tells the code which independent variable we used

% constraint - the physical constraint (string) - there are four
%       different string options for a physical constraint:
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
constraint = profile_type;              % string that tells the code which physical constraint to use



re = create_droplet_profile2([r_top, r_bottom], z, indVar, constraint);     % microns - effective radius vector


% Set the droplet distribution type
dist_str = GN_inputs.RT.drop_distribution_str;                                 % droplet distribution

% define the droplet distribution variance
% This should be the same length as re
% A distribution variance must be defined for each re value

% -- For now, lets assume this is constant --
dist_var = linspace(GN_inputs.RT.drop_distribution_var,GN_inputs.RT.drop_distribution_var, GN_inputs.RT.cloud_layers);              % distribution variance

vert_homogeneous_str = GN_inputs.RT.vert_homogeneous_str;     % This tells the function to create a multi-layered cloud
% define the boundaries of the cloud in Z-space
z_topBottom = [z(end), z(1)];                    % km - boundaries of the altitude vector. 

% Tell the code to use a pre-computed mie table for the extinction
% efficiency, or to use the value of the extinction paradox -> Qe = 2
parameterization_str = GN_inputs.RT.parameterization_str;


% -----------------------------------
% ---- Write a Water Cloud file! ----
% -----------------------------------

% ------------------------------------------------------
% --------------------VERY IMPORTANT ------------------
% ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
% ------------------------------------------------------
loop_var = 0;

wc_filename = write_wc_file(re,tau_c,z_topBottom, wavelength_tau_c(1,1), dist_str,...
    dist_var, vert_homogeneous_str, parameterization_str, loop_var);


% ------------------------------------------------------
% ------------------------------------------------------


% ----- Write an INP file --------
GN_names.inp = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename);
    
% now lets write the output names
    
GN_names.out = writeOutputNames(GN_names.inp);

% ---- Run uvspec for the files created -----
[measurement_estimate,~] = runReflectanceFunction_4modis_gaussNewton(GN_names,INP_folderName,saveCalculations_fileName, GN_inputs.spec_response);






end