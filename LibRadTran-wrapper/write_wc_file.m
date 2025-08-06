%% This function will write a .DAT water cloud file for LibRadTran

% This function will read and interpolate precomputed Mie calculations for
% water droplets of varrying radii.

% INPUTS:
%   (1) re - effective droplet radius (microns) - this is either a single
%   value, a vector, or a matrix. A single value for re tells the function
%   to create a cloud with a single layer containing a constant droplet
%   radius value. A vector tells the function to create a single wc file
%   with a droplet profile. The length of the vector is equal to the number
%   of layers modeled. A matrix tells the function to create multiple wc
%   files, where the number of columns is equal to the number of wc files
%   created. The number of rows is equal to the number of layers modeled.
%   To create many water cloud files at once that model a homogenous cloud,
%   simply set the column vectors of re to be identical values.
%   ***IMPORTANT*** the re values must start at the cloud bottom with the
%   first value (or row). The last value (or row) is the droplet size at
%   cloud top.

%   (2) tau_c - cloud optical depth (unitless) - this is the cloud optical
%   depth, which is a monochromatic calculation. There is a single value
%   that defines the optical depth of the cloud. LibRadTran defines cloud
%   files in terms of two values that do not depend on wavelength: re and
%   liquid water content (LWC). But usually people talk about clouds as
%   having a certain droplet size and a certain optical thickness. Enter a
%   single value for a single wc file, or a vector if you're creating
%   multiple wc files. If re is a matrix, and tau_c is a single value, each
%   wc file will have the entered tau_c. If each value in the re matrix
%   needs a unique tau value

%   (3) z_topBottom - altitude above sea level (kilometers) - this is a
%   vector with two values: [z_cloudTop, z_cloudBottom]. LibRadTran
%   constructs a cloud by creating layers, where each layer is homogenous
%   untill the next layer is defined. z_cloudTop defines where the cloud
%   ends; this is where the LWC should go to zero. This function will
%   compute a z vector equal in length to that of re using z_topBottom and
%   the geometric thickness H. If re is a matrix, the function expects
%   z_topBottom to be a matrix, where each column is a new wc file. If re
%   is a matrix, and z_topBottom is a vector, then this will be used for
%   every wc file.


%   (5) lambda - wavelength that defines the cloud optical depth
%   (nanometers) - If creating a single wc file, lambda is a single value. If
%   creating multiple wc files, lambda is a vector equal in length to the
%   number of columns of re. If re is a matrix and lambda is a single
%   value, this value will be used for each wc file created.

%   (6) distribution_str - a string telling the code which droplet size
%   distribution to use  - One can chose from two options:
%       (a) 'mono' - monodispersed distribution
%       (b) 'gamma' - gamma droplet distribution.
%       *** IMPORTANT *** For now, this function will NOT
%       use precomputed mie calculations using a gamma droplet
%       distribution. The values returned by LibRadTran appear erroneously
%       high. Instead, if one wishes to use a gamma droplet distribution,
%       the homogenous pre-computed mie table will be used to retrieve mie
%       properties, and then this code will integrate those values over the
%       size distribution.

%   (6) distribution_var - the variance of the size distribution, if
%   applicable. If one is modelling a homogenous cloud, this input will be
%   ignored, and the value can be anything.

%   (7) vert_homogeneity_str - a string telling the code if the cloud is to be
%   modeled as vertically homogeneous. If vertically homogenous, the code will
%   assume every single r_e value represents a single cloud with a constant
%   droplet radius. If vertically non-homogenous, each column of re is assumed
%   to be a single cloud with a droplet profile.
%       (a) 'vert-homogeneous' - homogenous cloud where the entire cloud can be
%       modeled as a single layer with constant properties. If this option
%       is chosen, the code will expect a vector for re, where each entry
%       represents a different cloud.
%       (b) 'vert-non-homogeneous' - a non-homogeneous cloud implies a cloud
%       with multiple layers where the properties vary. If this option is
%       chosen, the code expects a single column vector for re, or a
%       matrix, where each column vector represents a cloud


%   (8) parameterization_str - a string telling the code which
%   parameterization to use when using optical depth and droplet radius to
%   compute the liquid water content. There are 2 options
%       (a) 'mie' - this option uses a pre-computed mie table and
%       interpolates to find values for the extinction efficiency at each
%       wavelength and each droplet radius
%       (b) '2Limit' - this option assumes we are close to the extinction
%       paradox limit and sets the extinction efficiency at a constant
%       value of 2

%   (9) ind_var - the independent variable used to define the effective
%  radius profile. Options are
%           (1) 'optical_depth'
%           (2) 'altitude'

%   (10) compute_weighting_functions - a true or false flag that tells the
%   compute to compute weighting functions by creating N wc files where N
%   is equal to the number of cloud layers. The first file has all N
%   layers, and each successive file removes one cloud layer from the
%   bottom

%   (11) computer_name - the computer this code is running on

%   (12) index - this is the unique identifier that ensures files are not
%   written over one another. If one file is created, the index will be 1.
%   If many fiels are written in a for loop, each file will be tagged with
%   the number in the loop.

%   (13) num_re_parameters - There can be 1 free droplet size parameter for 
%   homogeneous clouds, or 2/3 for inhomogeneous clouds. For inhomogeneous,
%   there are either 2 free parameters (re_top and re_bot)
%   or 3 (re_top, re_middle, and re_bot)


% OUTPUTS:
%   (1) .Dat file saved in the libRadTran folder:
%   /.../libRadtran-2.0.4/data/wc

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

%%

function [fileName, lwc, ext_bulk_coeff_per_LWC] = write_wc_file(re, tau_c, z_topBottom, lambda, distribution_type,...
    distribution_var, vert_homogeneous_str, parameterization_str, ind_var, compute_weighting_functions,...
    computer_name, index, num_re_parameters)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 8 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=13
    error([newline,'Should be 13 inputs: droplet effective radius, optical depth, altitude,',...
        [' wavelength, droplet distribution type, variance of the droplet distribution,',...
        ' homogeneity type, the parameterization used to compute LWC, the independent vertical variable',...
        ', a flag telling the code to compute weighting functions, the computer name,',...
        'a unique file index, and a number associated with the re free parameters.'], newline])
end

% Check to make sure re is the same length as the altitude vector
[r_re,c_re] = size(re);

% first check to see if z_topBottom is a vector or a matrix
if size(z_topBottom,1)==1 || size(z_topBottom,2)==1
    % If true, then there must be two entries
    if length(z_topBottom)~=2
        error([newline,'Need two values for z_topBottom: altitude at cloud bottom top and cloud bottom', newline])
    end

    % make sure its a  column vector
    z_topBottom = reshape(z_topBottom,[],1);
    H = z_topBottom(1) - z_topBottom(2);            % geometric thickness

elseif size(z_topBottom,1)>1 && size(z_topBottom,2)>1
    % if true, then there can only be two rows, and it must be equal in
    % size to the r matrix
    if size(z_topBottom,1)~=2
        error([newline,'Need two values for z_topBottom: altitude at cloud bottom top and cloud bottom', newline])

    elseif size(z_topBottom,2)~=size(re,2) || size(z_topBottom,2)==1
        error([newline,'z_topBottom must have the same number of columns as re, or a single column that is used for all wc files.', newline])

    end

    H = z_topBottom(:,1) - z_topBottom(:,2);        % geometric cloud thickness

end



if length(lambda)>1 && length(lambda)~=size(re,2)

    error([newline,'Lambda must be either a single value or a vector equal in legnth to the number of columns in re.', newline])
end

if length(tau_c)>1 && length(tau_c)~=size(re,2)

    error([newline,'The optical depth must be either a single value or a vector equal in legnth to the number of columns in re.', newline])
end


% Check to make sure the distribution string is one of two possible values

if strcmp(distribution_type, 'mono')==false && strcmp(distribution_type, 'gamma')==false

    error([newline,'I dont recognize the droplet distribution. Must be either "mono" or "gamma"', newline])
end


% Check to make sure the homogeneity string is one of two possible values

if strcmp(vert_homogeneous_str, 'vert-homogeneous')==false && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==false

    error([newline,'I dont recognize the homogeneity string. Must be either "vert-homogenous" or "vert-non-homogeneous"', newline])
end


% ----- Check to see if there are any NaNs in the radius vector -----

if any(isnan(re))==true

    error([newline, 'The effective radius has atleast one NaN value.', newline])
end







% Lets set up a few warnings incase the values of effective radius are
% outside the bounds of the Hu and Stamnes or Mie Interpolate
% parameterization
if size(re,1)>1 && size(re,2)>1
    if strcmp(parameterization_str, 'hu')==true && any(any(re<2.5)) || any(any(re>60))
        warning([newline, 'At least one value in r_{e} is outside the range of the Hu and Stamnes parameterization',newline,...
            'This is the default parameterization used to convert water cloud parameters to optical properites.',newline,...
            'The range of acceptable values for this parameterization is [2.5, 60] microns.',newline]);
    end

    if strcmp(parameterization_str, 'mie')==true && any(any(re<1)) || any(any(re>25))
        warning([newline,'At least one value in r_{e} is greater than 25 microns, which is the upper limit',...
            ' to the Mie Interpolate parameterization that computes optical properties from the water cloud ',...
            'parameters. The netcdf file downloaded from LibRadTrans website includes values of re up to',...
            ' 25 microns. The acceptable range for Mie Interpolation is [1, 25] microns.',newline]);
    end

else

    if strcmp(parameterization_str, 'hu')==true && any(re<2.5) || any(re>60)
        warning([newline, 'At least one value in r_{e} is outside the range of the Hu and Stamnes parameterization',newline,...
            'This is the default parameterization used to convert water cloud parameters to optical properites.',newline,...
            'The range of acceptable values for this parameterization is [2.5, 60] microns.',newline]);
    end

    if strcmp(parameterization_str, 'mie')==true && any(any(re<1)) || any(any(re>25))
        warning([newline,'At least one value in r_{e} is greater than 25 microns, which is the upper limit',...
            ' to the Mie Interpolate parameterization that computes optical properties from the water cloud ',...
            'parameters. The netcdf file downloaded from LibRadTrans website includes values of re up to',...
            ' 25 microns. The acceptable range for Mie Interpolation is [1, 25] microns.',newline]);
    end

end






% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    water_cloud_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/';

elseif strcmp(computer_name,'andrewbuggee')==true

    water_cloud_folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/';

elseif strcmp(computer_name,'curc')==true

    water_cloud_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/wc/';


end



%%

% ------------------------------------------------------------
% ---------------------- COMPUTE LWC -------------------------
% ------------------------------------------------------------


%rho_liquid_water = 997048;                 % grams/m^3 - density of liquid water at 0 C - Wolfram Alpha
rho_liquid_water = 1e6;                     % grams/m^3
% --------------------------------------------------


% --- STEP THROUGH COLUMNS OF re -----

% if re is a matrix, the code assumes that each column vector is a single
% cloud with a vertically-inhomogenous droplet profile defined by the
% values in the column
if size(re,1)>1 && size(re,2)>1

    % We want to get all the mie properties we need before we loop through
    % the number of files

    num_files_2write = size(re,2);          % number of wc files to create



    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing water cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values



    if strcmp(parameterization_str, 'mie')==true

        % lambda has to be the same size as re
        if length(lambda)==num_files_2write

            % remember that vectorizing a matrix (:) will stack column vectors
            % on top of one another

            lambda = reshape(lambda,1,[]); % each column represents a new file

            % create identical column vecotrs
            lambda = repmat(lambda,size(re,1),1);


            % **** ONLY INTERPOLATING HOMOGENOUS MIE COMPUTATIONS ****
            % ********************************************************
            % IF GAMMA DISTRIBUTION DESIRED, CODE WILL MANUALLY
            % INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED

            if strcmp(distribution_type,'gamma')==true

                % If we wish to estimate the mie properties of liquid water
                % for a distribution of droplets, then we can skip the
                % pre-computed mie tables and estimate the values using the
                % average_mie_over_size_distribution directly

                % This function only deals with liquid water clouds
                % define the index of refraction
                index_of_refraction = 'water';

                % integrate over a size distribution to get an average
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
                    index_of_refraction, distribution_type, computer_name, index);

            elseif strcmp(distribution_type,'mono')==true
                yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

            else

                error([newline,'Invaled distribution type',newline])

            end



        elseif isscalar(lambda)

            % **** ONLY INTERPOLATING HOMOGENOUS MIE COMPUTATIONS ****
            % ********************************************************
            % IF GAMMA DISTRIBUTION DESIRED, CODE WILL MANUALLY
            % INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED

            if strcmp(distribution_type,'gamma')==true

                % If we wish to estimate the mie properties of liquid water
                % for a distribution of droplets, then we can skip the
                % pre-computed mie tables and estimate the values using the
                % average_mie_over_size_distribution directly

                % This function only deals with liquid water clouds
                % define the index of refraction
                index_of_refraction = 'water';

                % integrate over a size distribution to get an average
                [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
                    index_of_refraction, distribution_type, computer_name, index);

            elseif strcmp(distribution_type,'mono')==true
                yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

            else

                error([newline,'Invaled distribution type',newline])

            end

        end

    elseif strcmp(parameterization_str,'2limit')==true
        % set the value to be 2 for all calculations
        yq = 2*ones(length(re),5);

    end






    % if re is a vector and the vertically homogeneity is defined as
    % 'vert-non-homogeneous' then the the code assumes the vector defines a
    % single cloud with a droplet profile
elseif ((size(re,1)==1 && size(re,2)>1) || (size(re,1)>1 && size(re,2)==1)) &&...
        strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true

    num_files_2write = 1;

    % re must be a column vector
    re = reshape(re,[],1);

    % if the independent variable is optical depth, we need to flip this
    % vector. This function works in z space, where z increases from bottom
    % to top, therefore we need to start with the radius at max optical
    % depth
    if strcmp(ind_var, 'optical_depth')==true
        re = flipud(re);
    end

    % the distribution variance should be a column vector
    distribution_var = reshape(distribution_var, [], 1);

    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing water cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values

    if strcmp(parameterization_str,'mie')==true


        if strcmp(distribution_type,'gamma')==true

            % -------------------------------------------------------------
            % --- MANUALLY INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED ---
            % -------------------------------------------------------------
            % If we wish to estimate the mie properties of liquid water
            % for a distribution of droplets, then we can skip the
            % pre-computed mie tables and estimate the values using the
            % average_mie_over_size_distribution directly

            % This function only deals with liquid water clouds
            % define the index of refraction
            %             index_of_refraction = 'water';

            % this loop applies to a vertical droplet profile. For now we
            % will apply the same distribution variance to each level in
            % the cloud.

            % integrate over a size distribution to get an average
            %             [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var,...
            %                 lambda,index_of_refraction, distribution_type, index);


            % -------------------------------------------------------
            % ----------- USE LIBRADTRAN MIE CALCULATIONS -----------
            % -------------------------------------------------------
            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For water, since the density is 1 g/m^3, we
            % can simply multiply the output with the liquid water content
            % and integrate over the path to get the optical depth.


            % What mie code should we use to compute the scattering properties?
            mie_program = 'MIEV0';               % type of mie algorithm to run

            % This function only deals with liquid water clouds
            % define the index of refraction
            index_of_refraction = 'water';

            size_distribution = {'gamma', distribution_var(1)};           % droplet distribution

            % Do you want a long or short error file?
            err_msg_str = 'verbose';



            % The radius input is defined as [r_start, r_end, r_step].
            % where r_step is the interval between radii values (used only for
            % vectors of radii). A 0 tells the code there is no step. Finally, the
            % radius values have to be in increasing order.
            ext_bulk_coeff_per_LWC = zeros(length(re), 1);

            for rr = 1:length(re)

                mie_radius = [re(rr), re(rr), 0];    % microns


                % Create a mie file
                [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
                    mie_radius, lambda, size_distribution, err_msg_str, computer_name, index, round(re(rr), 4));

                % run the mie file
                [~] = runMIE(mie_folder,input_filename,output_filename, computer_name);

                % Read the output of the mie file
                [ds,~,~] = readMIE(mie_folder,output_filename);

                ext_bulk_coeff_per_LWC(rr) = ds.Qext;       % km^-1 / (cm^3 / m^3)

            end
            % --------------------------------------------------------------




        elseif strcmp(distribution_type,'mono')==true
            yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

        else

            error([newline,'Invaled distribution type',newline])

        end

    elseif strcmp(parameterization_str, '2limit')==true

        yq = 2*ones(length(re),5);

    end




    % When creating weighting functions, we want to build a cloud layer by
    % layer. So it may be vertically inhomogeneous but the first layer is not,
    % because by definition a single layer cloud is homogeneous.
elseif isscalar(re) && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true


    % treat this single layer cloud as vertically homogeneous

    num_files_2write = 1;


    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing water cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values

    if strcmp(parameterization_str,'mie')==true




        if strcmp(distribution_type,'gamma')==true

            % -------------------------------------------------------------
            % --- MANUALLY INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED ---
            % -------------------------------------------------------------

            % If we wish to estimate the mie properties of liquid water
            % for a distribution of droplets, then we can skip the
            % pre-computed mie tables and estimate the values using the
            % average_mie_over_size_distribution directly

            % This function only deals with liquid water clouds
            % define the index of refraction
            %             index_of_refraction = 'water';
            %
            %             % integrate over a size distribution to get an average
            %             [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
            %                 index_of_refraction, distribution_type, index);
            % -------------------------------------------------------------


            % -------------------------------------------------------
            % ----------- USE LIBRADTRAN MIE CALCULATIONS -----------
            % -------------------------------------------------------
            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For water, since the density is 1 g/m^3, we
            % can somply multiply the output with the liquid water content
            % and integrate over the path to get the optical depth.


            % What mie code should we use to compute the scattering properties?
            mie_program = 'MIEV0';               % type of mie algorithm to run

            % This function only deals with liquid water clouds
            % define the index of refraction
            index_of_refraction = 'water';

            size_distribution = {'gamma', distribution_var(1)};           % droplet distribution

            % Do you want a long or short error file?
            err_msg_str = 'verbose';


            % The radius input is defined as [r_start, r_end, r_step].
            % where r_step is the interval between radii values (used only for
            % vectors of radii). A 0 tells the code there is no step. Finally, the
            % radius values have to be in increasing order.

            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For water, since the density is 1 g/m^3, we
            % can somply multiply the output with the liquid water content
            % and integrate over the path to get the optical depth.

            mie_radius = [re, re, 0];    % microns


            % Create a mie file
            [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
                mie_radius, lambda, size_distribution, err_msg_str, computer_name, index, round(re, 4));

            % run the mie file
            [~] = runMIE(mie_folder,input_filename,output_filename, computer_name);

            % Read the output of the mie file
            [ds,~,~] = readMIE(mie_folder,output_filename);

            ext_bulk_coeff_per_LWC = ds.Qext;       % km^-1 / (cm^3 / m^3)

            % --------------------------------------------------------------





        elseif strcmp(distribution_type,'mono')==true
            yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

        else

            error([newline,'Invaled distribution type',newline])

        end


    elseif strcmp(parameterization_str, '2limit')==true

        yq = 2*ones(length(re),5);

    end






    % if the re input is a vector and the homogenous string is defined as
    % vertically homogeneous, then the code assumes each value in the vector is
    % a single cloud, and the each value defines the homogenous droplet size
    % for that cloud.
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true

    num_files_2write = 1;


    % -------------------------------------------------------------------
    % ------ open the precomputed mie table and interpolate! ------------
    % -------------------------------------------------------------------

    % for writing water cloud files, we only need the extinction efficiency
    % Since this function is used often, we've created a file with just Q_ext

    justQ = true;                       % Load the precomputed mie table of only Q_ext values

    if strcmp(parameterization_str,'mie')==true




        if strcmp(distribution_type,'gamma')==true

            % -------------------------------------------------------------
            % --- MANUALLY INTEGRATE OVER THE SIZE DISTRIBUTION DEFINED ---
            % -------------------------------------------------------------

            % If we wish to estimate the mie properties of liquid water
            % for a distribution of droplets, then we can skip the
            % pre-computed mie tables and estimate the values using the
            % average_mie_over_size_distribution directly

            % This function only deals with liquid water clouds
            % define the index of refraction
            %             index_of_refraction = 'water';
            %
            %             % integrate over a size distribution to get an average
            %             [~, Qe_avg, ~] = average_mie_over_size_distribution(re, distribution_var, lambda,...
            %                 index_of_refraction, distribution_type, index);
            % -------------------------------------------------------------


            % -------------------------------------------------------
            % ----------- USE LIBRADTRAN MIE CALCULATIONS -----------
            % -------------------------------------------------------
            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For water, since the density is 1 g/m^3, we
            % can somply multiply the output with the liquid water content
            % and integrate over the path to get the optical depth.


            % What mie code should we use to compute the scattering properties?
            mie_program = 'MIEV0';               % type of mie algorithm to run

            % This function only deals with liquid water clouds
            % define the index of refraction
            index_of_refraction = 'water';

            size_distribution = {'gamma', distribution_var(1)};           % droplet distribution

            % Do you want a long or short error file?
            err_msg_str = 'verbose';


            % The radius input is defined as [r_start, r_end, r_step].
            % where r_step is the interval between radii values (used only for
            % vectors of radii). A 0 tells the code there is no step. Finally, the
            % radius values have to be in increasing order.

            % Libradtran doesn't compute the efficieny when a distribution
            % is specified. It computes the bulk coefficient per unit
            % concentration. For water, since the density is 1 g/m^3, we
            % can somply multiply the output with the liquid water content
            % and integrate over the path to get the optical depth.
            ext_bulk_coeff_per_LWC = zeros(length(re), 1);

            for rr = 1:length(re)

                mie_radius = [re(rr), re(rr), 0];    % microns


                % Create a mie file
                [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
                    mie_radius, lambda, size_distribution, err_msg_str, computer_name, index, round(re(rr), 4));

                % run the mie file
                [~] = runMIE(mie_folder,input_filename,output_filename, computer_name);

                % Read the output of the mie file
                [ds,~,~] = readMIE(mie_folder,output_filename);

                ext_bulk_coeff_per_LWC(rr) = ds.Qext;       % km^-1 / (cm^3 / m^3)

            end
            % --------------------------------------------------------------





        elseif strcmp(distribution_type,'mono')==true
            yq = interp_mie_computed_tables([repmat(lambda,numel(re),1), re], 'mono', justQ);

        else

            error([newline,'Invaled distribution type',newline])

        end


    elseif strcmp(parameterization_str, '2limit')==true

        yq = 2*ones(length(re),5);

    end





else

    error([newline,'re is not a vector or a matrix! Check your inputs!', newline])

end



% grab the extinction efficiency values

if strcmp(distribution_type,'gamma')==true
    %     Qext = Qe_avg';         % Extinction efficiency
    %Qext = linspace(2.0816, 2.0816, length(re))';        % value to match libRadTran

elseif strcmp(distribution_type,'mono')==true
    Qext = reshape(yq(:,3),[],num_files_2write);         % convert this back into a matrix corresponging to re

end




% now we will step through each wc file that needs to be created

% if H is a single value and num_files_2write is greater than 1, we will
% repeat it to create a vector with the same length
if isscalar(H) && num_files_2write>1
    H = repmat(H,num_files_2write,1);           % km - geometric thickness
end


% if tau_c is a single value and num_files_2write is greater than 1, we will
% repeat it to create a vector with the same length
if isscalar(tau_c) && num_files_2write>1
    tau_c = repmat(tau_c,num_files_2write,1);
end

if size(z_topBottom,2)==1 && num_files_2write>1
    z_topBottom = repmat(z_topBottom,1,num_files_2write);
end

% How many files are being created?
fileName = cell(1,num_files_2write);


% How many layers to model in the cloud?
if size(re,1)>1 && size(re,2)>1 && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true
    nLayers = size(re,1)+1;             % Number of altitude levels we need to define a cloud
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true
    nLayers = length(re)+1;             % Number of altitude levels we need to define a cloud
elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true
    nLayers = 1;

end


for nn = 1:num_files_2write


    % -------------------------------------------
    % ------ Create altitude vector! ------------
    % -------------------------------------------

    % the length of the altitude vector should be 1 unit longer than the length
    % of the effective radius. Thats because the last value in the altitude
    % vector is the altitude at cloud top, where the LWC has gone to zero


    % z must be a column vector
    if nLayers==1
        z = flipud(z_topBottom(:,nn));

    else
        z = linspace(z_topBottom(2,nn), z_topBottom(1,nn), nLayers)';                 % km - altitude vector
    end


    % -------------------------------------------------------------------
    % ------ compute number concentration and liquid water content ------
    % -------------------------------------------------------------------

    % we need a number concentration for each file that is created

    if nLayers>1 && strcmp(distribution_type, 'gamma')==true

        % We could just integrate the size distributiion to get the total
        % number concentration, but we've chosen an arbitrary value for the
        % effective variance because it doesn't have much effect on
        % reflectance measurements over the solar spectral region. More
        % importantly, we want to connect two user defined variables, cloud
        % optical depth and effective radius, to the number concentration,
        % and thus the liquid water content.
        %         z_meters = (z(1:end-1)-z(1))*1e3;       % meters - geometric depth, normalized
        %         re_meters = (re(:,nn)*1e-6);            % meters - effective radius converted to meters
        %
        %         Nc = tau_c(nn)./(pi*trapz(z_meters, Qext(:,nn).* re_meters.^2));                % m^(-3) - number concentration
        %
        %         % ------------------------------------------------------------------
        %         % --- Solve for the total Liquid Water Content over the entire cloud ---
        %         % number concentration is constant with height. We make the
        %         % assumption that all droplets can be modeled as the effective
        %         % radius. So the LWC simple changes with effective radius
        %         % ** LibRadTran requires LWC in units of grams/m^3 **
        %         lwc = 4/3 * pi * rho_liquid_water * re_meters.^3 .* Nc;                    % g/m^3 - grams of water per meter cubed of air
        % -----------------------------------------------------------------



        % ----------------------------------------------------------------
        % ******** Integrating over monodispersed mie caluclation ********
        % ----------------------------------------------------------------
        % -- Assuming liquid water content increases linearly with depth -
        %         re_meters = (re(:,nn)*1e-6);            % meters - effective radius converted to meters
        %         z_meters_midpoint = ((z(1:end-1)-z(1)) + (z(2)-z(1))/2)*1e3;       % meters - geometric depth, normalized
        %         dz = z_meters_midpoint(2)-z_meters_midpoint(1);           % meters
        %
        %
        %         slope = (4*rho_liquid_water * tau_c) /(3*dz * sum(Qext .* z_meters_midpoint ./re_meters));     % g/m^3/m - slope of the lwc profile
        %
        %         % solve for the linear liquid water content profile
        %         lwc = slope * z_meters_midpoint;                     % g/m^3 - grams of water per meter cubed of air
        % ----------------------------------------------------------------



        % -----------------------------------------------------------------
        % ** Using libRadTran mie calculations with a size distribution ***
        % -----------------------------------------------------------------
        % ** Assuming liquid water content increases linearly with depth **

        %z_kilometers_midpoint = ((z(1:end-1)-z(1)) + (z(2)-z(1))/2);       % kilometers - geometric depth at midpoint of each layer
        z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
        dz_km = z(2) - z(1);           % kilometers

        %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
        slope = tau_c(nn) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/m - slope of the lwc profile

        % solve for the linear liquid water content profile
        %lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
        lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air
        % ----------------------------------------------------------------



        % -----------------------------------------------------------------
        % ******** compute LWC by integrated the size distribution ********
        % -----------------------------------------------------------------
        % *** There is another way to solve for the LWC ***
        % We've made an assumption about the droplet size distribution and
        % we've computed the total number concentration. We can solve for
        % the liquid water content by integrating the size distribution
        % *** IMPORTANT *** We have to play with the distribution width to
        % get the correct optical depth
        %         if strcmp(distribution_str,'gamma')==true
        %
        %             %distribution_var = 27;
        %             lwc = zeros(size(re));
        %
        %             for zz = 1:length(re)
        %
        %                 [nr,r] = gamma_size_distribution_libRadTran2(re(zz), distribution_var(zz), Nc);       % [#/micron/m^3 , microns] - gamma droplet size distribution
        %                 lwc(zz) = trapz( r , 4/3 * pi * rho_liquid_water * (r*1e-6).^3 .* nr);                % g/m^3 - grams of water per meter cubed of air
        %
        %             end
        %
        %         end
        % ------------------------------------------------------------------
        % ------------------------------------------------------------------


        % create the water cloud file name
        if index==0

            if num_re_parameters==2

                % print top and bottom

                fileName{nn} = ['WC_rtop',num2str(round(re(end,nn),3)),'_rbot',num2str(round(re(1,nn),3)),'_T',num2str(round(tau_c(nn),3)),...
                    '_', distribution_type,'_nn',num2str(nn), '.DAT'];

            elseif num_re_parameters==3

                % print top, middle and bottom
                % find the cloud middle
                [~, idx_middle] = min(abs(z - ((max(z) - min(z))/2 + min(z)) ) );

                fileName{nn} = ['WC_rtop',num2str(round(re(end,nn),3)), '_rmid', num2str(round(re(idx_middle,nn),3)),...
                    '_rbot',num2str(round(re(1,nn),3)),'_T',num2str(round(tau_c(nn),3)),...
                    '_', distribution_type,'_nn',num2str(nn), '.DAT'];

            end


        elseif index>0

            if num_re_parameters==2

                % print top and bottom

            fileName{nn} = ['WC_rtop',num2str(round(re(end,nn),3)),'_rbot',num2str(round(re(1,nn),3)),'_T',num2str(round(tau_c(nn),3)),...
                '_', distribution_type,'_nn',num2str(index), '.DAT'];

            elseif num_re_parameters==3

                % print top, middle and bottom
                % find the cloud middle
                [~, idx_middle] = min(abs(z - ((max(z) - min(z))/2 + min(z)) ) );

                fileName{nn} = ['WC_rtop',num2str(round(re(end,nn),3)), '_rmid',num2str(round(re(idx_middle,nn),3)),...
                    '_rbot',num2str(round(re(1,nn),3)),'_T',num2str(round(tau_c(nn),3)),...
                '_', distribution_type,'_nn',num2str(index), '.DAT'];

            end

        end



    elseif nLayers==1 && strcmp(distribution_type, 'gamma')==true

        % ----------------------------------------------------------------
        % --- If there is one cloud layer, this is a homogensous cloud ---
        % ----------------------------------------------------------------

        %Nc = tau_c(nn)./(pi*(H(nn)*1e3)*Qext(nn).*(re(nn)*1e-6).^2);                 % m^(-3) - number concentration

        % Compute Liquid Water Content
        %lwc = 4/3 * pi * rho_liquid_water * (re(nn)*1e-6).^3 .* Nc;                  % g/m^3 - grams of water per meter cubed of air

        % Compute Liquid Water Content
        lwc = tau_c(nn)./(ext_bulk_coeff_per_LWC(nn) .* H(nn));                           % g/m^3 - grams of water per meter cubed of air


        % create the water cloud file name
        if index==0
            fileName{nn} = ['WC_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(nn), '.DAT'];
        elseif index>0

            fileName{nn} = ['WC_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(index), '.DAT'];
        end




    elseif nLayers==1 && strcmp(distribution_type, 'mono')==true

        % ----------------------------------------------------------------
        % --- If there is one cloud layer, this is a homogensous cloud ---
        % ----------------------------------------------------------------

        %Nc = tau_c(nn)./(pi*(H(nn)*1e3)*Qext(nn).*(re(nn)*1e-6).^2);                 % m^(-3) - number concentration

        % Compute Liquid Water Content
        %lwc = 4/3 * pi * rho_liquid_water * (re(nn)*1e-6).^3 .* Nc;                  % g/m^3 - grams of water per meter cubed of air

        % Compute Liquid Water Content
        lwc = 4/3 * (re(nn)*1e-6) * rho_liquid_water * tau_c(nn)./...
            (Qext(nn) * (H(nn)*1e3));                                                % g/m^3 - grams of water per meter cubed of air

        % create the water cloud file name
        if index==0
            fileName{nn} = ['WC_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(nn), '.DAT'];
        elseif index>0

            fileName{nn} = ['WC_r',num2str(round(re(nn))),'_T',num2str(round(tau_c(nn))),'_', distribution_type,...
                '_nn',num2str(index), '.DAT'];
        end





    end




    % ------------------------------------------------------------
    % ----------------- WE NEED TO APPEND ZEROS ------------------
    % ------------------------------------------------------------

    % Wherever the cloud is, there needs to be zeros at the cloud top altitude,
    % and below the cloud bottom altitude. This information tells LibRadTran
    % where the boundaries of the cloud are

    % both the effective radius and the LWC need zeros on either boundary,
    % unless if the cloud is at the surface

    if (r_re==1 || c_re==1) && strcmp(vert_homogeneous_str, 'vert-non-homogeneous')==true

        if z(1)==0
            % If true, then the cloud starts at the surface and we only append
            % zeros above the cloud
            re_2write = [re; 0];
            lwc_2write = [lwc; 0];
            z_2write = z;

        else
            % In this case, we need zeros below the cloud bottom, and at cloud
            % top
            z_2write = [0; z];                 % create a value at the surface where the cloud parameters go to zero
            re_2write = [0; re; 0];
            lwc_2write = [0; lwc; 0];

        end

    elseif (size(re,1)==1 || size(re,2)==1) && strcmp(vert_homogeneous_str, 'vert-homogeneous')==true

        if z(1)==0
            % If true, then the cloud starts at the surface and we only append
            % zeros above the cloud
            re_2write = [re(nn); 0];
            lwc_2write = [lwc; 0];
            z_2write = z;

        else
            % In this case, we need zeros below the cloud bottom, and at cloud
            % top
            z_2write = [0; z];                 % create a value at the surface where the cloud parameters go to zero
            re_2write = [0; re(nn); 0];
            lwc_2write = [0; lwc; 0];

        end

    elseif (size(re,1)>1 && size(re,2)>1)

        % Cloud top height defines the altitude where there is no cloud.

        % if the minimum z value is 0 then the cloud is at the surface
        if z(1)==0
            % then we only append zeros above the cloud
            re_2write = [re(:,nn); 0];
            lwc_2write = [lwc; 0];
            z_2write = z;

        else
            % Then we need zeros on either end
            z_2write = [0; z];
            re_2write = [0; re(:,nn); 0];
            lwc_2write = [0; lwc; 0];

        end

    end




    % ------------------------------------------------------------
    % ---------------------- WRITE WC FILE -----------------------
    % ------------------------------------------------------------


    if compute_weighting_functions==true

        % to compute the weighting functions using the methods of Platnick
        % (2000), we need to compute the reflectance at cloud top (or TOA)
        % by incrementally adding a layer at cloud bottom. Using the z,re,
        % and lwc above, we will iterate through each layer and create n
        % files where n is equal to the number of cloud layers. Each
        % iteration will remove a layer from cloud bottom.

        % compute the optical depth at each layer
        %tau_layer = lwc .* ext_bulk_coeff_per_LWC .* dz_km;       % optical thickness of each homogeneous layer
        % we need a filename for each layer
        filename = cell(nLayers-1, 1);

        for LL = 1:nLayers-1

            % define the filename
            filename{LL} = [fileName{1}(1:end-7), 'layers1-', num2str(nLayers-LL),...
                '.DAT'];
            % Create the water cloud file
            fileID = fopen([water_cloud_folder_path,filename{LL}], 'w');

            % fprintf writes lines in our text file from top to botom
            % wc.DAT files are written with the higher altitudes at the top, and the
            % surface at the bottom

            % to write column vectors in a text file, we have to store them as row
            % vectors

            toWrite = [flipud(z_2write([1, LL+1:end]))'; flipud(lwc_2write([1, LL+1:end]))';...
                flipud(re_2write([1, LL+1:end]))'];

            % Create the opening comment lines of the WC.DAT file

            fprintf(fileID, '%s %10s %7s %8s \n','#','z','LWC','R_eff');
            fprintf(fileID, '%s %10s %7s %8s \n','#','(km)','(g/m^3)','(micron)');

            % Write in the data
            fprintf(fileID,'%12.3f %7.4f %8.3f \n', toWrite);
            fclose(fileID);


        end

        % redefine the output file names
        fileName = filename;


    else



        % Create the water cloud file
        fileID = fopen([water_cloud_folder_path,fileName{nn}], 'w');

        % fprintf writes lines in our text file from top to botom
        % wc.DAT files are written with the higher altitudes at the top, and the
        % surface at the bottom

        % to write column vectors in a text file, we have to store them as row
        % vectors

        toWrite = [flipud(z_2write)'; flipud(lwc_2write)'; flipud(re_2write)'];

        % Create the opening comment lines of the WC.DAT file

        fprintf(fileID, '%s %10s %7s %8s \n','#','z','LWC','R_eff');
        fprintf(fileID, '%s %10s %7s %8s \n','#','(km)','(g/m^3)','(micron)');

        % Write in the data
        fprintf(fileID,'%12.3f %7.4f %8.3f \n', toWrite);
        fclose(fileID);


    end

end


end
