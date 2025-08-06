%% This function will read a solar_flux file from LibRadTran

% INPUTS:
%   (1) wavelength - vector containing the wavelength range of interest
%   (nm) - user will provide a wavelength range in a vector with two
%   quantities [min, max]. The function will find all solar flux values
%   within this range.

%   (2) file_name - name of uvspec solar flux file (string) - there are
%   three options to choose from:
%       (a) 'kurudz_1.0nm.dat' - These data range from 250 - 10000 nm and
%       are spaced by 1 nm. These data were taken from LBLRTM 5.21
%       (http://www.meto.umd.edu/~bobe/LBLRTM/). The original Kurudz [1992]
%       data were converted to mW / (m2 nm) and averaged over 1nm intervals
%       centered around the given wavelength.

%       (b) 'kurudz_0.1nm.dat' - These data range from 250 - 10000 nm and
%       are spaced by 0.1 nm. These data were taken from LBLRTM 5.21
%       (http://www.meto.umd.edu/~bobe/LBLRTM/). The original Kurudz [1992]
%       data were converted to mW / (m2 nm) and averaged over 0.1nm intervals
%       centered around the given wavelength.

%       (c) 'atlas_plus_modtran' - Use this if you wish to have solar
%       flux values from 200 - 250 nm. These data range from 200 to 800 nm
%       with 0.05 nm resolution. You do not need an extension if using this
%       file!

%       (d)
%       'hybrid_reference_spectrum_p025nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.005 nm.  The original
%       data were converted to mW / (m2 nm)

%       (e)
%       'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.025 nm.  The original
%       data were converted to mW / (m2 nm)


%       (f)
%       'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat'
%       - this is a hybrind reference spectrum downloaded from LASP's
%       LISIRD tool (https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm)
%       These data range from 202 to 2730 nm
%       These data have a smapling resolution of 0.1 nm.  The original
%       data were converted to mW / (m2 nm)


% OUTPUTS:
%   (1) solar_flux (W/nm/m^2) - a vector containing the solar flux values
%   at each wavelength specified in the wavelength vector

%   (2) wavelength (nm) - a vector containing wavelengths where the solar flux is defined


% NOTES:
%   LibRadTran comes supplied with a folder of solar flux both measured and
%   modeled. This function allows the user to grab values from these solar
%   flux files using a specific wavelength range.

% By Andrew John Buggee
%%

function [solar_flux, wavelength] = read_solar_flux_file(wavelength_boundaries, file_name)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 2 inputs, wavelength bounds and a file_name
% to read


if nargin~=2
    error([newline,'Not enough inputs. Need 2: wavelength bounds and a filename', newline])
end

% Check to make sure wavelength has a length of 2

if length(wavelength_boundaries)==2

else
    error([newline,'The wavelength input is a vector consisting of: [wavelength_min, wavelength_max]', newline])
end


% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'andrewbuggee')==true


    solar_source_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/solar_flux/'];

elseif strcmp(computer_name,'anbu8374')==true

    %error('You havent stored the mie calculations on you laptop yet!')
    solar_source_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/';

elseif strcmp(computer_name, 'curc')==true

    % --- super computer ---
    solar_source_folder = '/projects/anbu8374/software/libRadtran-2.0.5/data/solar_flux/';

end


%%

% ------------------------------------------------------------
% ---------------------- Read Source File --------------------
% ------------------------------------------------------------



if strcmp(file_name,'kurudz_1.0nm.dat')

    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [250, 10000];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of kurudz_1.0nm.dat. Must be between [250, 10000] nm.', newline])
    end


    % ---- Open the File ----

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading

    % ------------------------------------------------------
    % -------- Reading .dat file using fscanf --------------
    % ------------------------------------------------------

    % define how the data is written
    % '%' starts a new character
    % '%*s' tells the code to skip string characters
    % 'f' tells us to look for floating point numbers
    % 'e' tells the code to look for exponential numbers
    % '5.1' tells us to look for a floating point number with 5 entities, 1
    % of which comes after the decimal point

    %     format_spec = '  %f %e';
    %     shape_output = [2 Inf];
    %
    %     % lets skip the first 11 lines, since they are commented
    %     for ii = 1:11
    %         fgets(file_id);
    %     end
    %
    %     % now the file pointer will be at the data
    %     A = fscanf(file_id, format_spec, shape_output)'; % sxtract data!

    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;



elseif strcmp(file_name,'kurudz_0.1nm.dat')

    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [250, 10000];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of kurudz_0.1nm.dat. Must be between [250, 10000] nm.', newline])
    end


    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;

elseif strcmp(file_name,'atlas_plus_modtran')

    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [200, 800];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of atlas_lus_modtran.txt. Must be between [200, 800] nm.', newline])
    end

    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;


elseif strcmp(file_name, 'hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc')==true

    % Use the TSIS-1 hypbring solar reference spectrum (W/m^2/nm)
    error([newline, 'No .dat file yet! Create one!', newline])

    % open the netCDF file
    %info = ncinfo([solar_source_folder, file_name]);

    % read the wavelength grid
    wavelength_grid = ncread([solar_source_folder, file_name], 'Vacuum Wavelength');     % nm

    % find the indices for the wavelengths that lie within out wavelengths
    % of interest
    index_wavelength = wavelength_grid>=wavelength_boundaries(1) & wavelength_grid<=wavelength_boundaries(2);

    % read in the solar spectral flux
    solar_flux_data = ncread([solar_source_folder, file_name], 'SSI');             % W/m^2/nm


    % store the values that we wish to keep
    wavelength = wavelength_grid(index_wavelength);
    solar_flux = solar_flux_data(index_wavelength);

    %     % Lets check to make sure the wavelength input is within bounds of the
    %     % file selected
    %
    %     wavelength_regime = [202, 2730];            % nanometers - wavelength boundaries
    %
    %     if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
    %         error([newline, 'Wavelength is out of the range of atlas_lus_modtran.txt. Must be between [200, 800] nm.', newline])
    %     end
    %
    %     % ------------------------------------------------------
    %     % -------- Reading .dat file using textscan ------------
    %     % ------------------------------------------------------
    %     % Or we could use the textscan() function instead, which allows us to define comments to ignore
    %
    %     file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading
    %
    %
    %     format_spec = '%f %f';                                  % two floating point numbers
    %     B = textscan(file_id, format_spec, 'CommentStyle','#');
    %
    %     index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);
    %
    %     wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range
    %
    %     solar_flux = B{2}(index_wavelength);                % Watts/m^2/nm - flux values at the corresponding wavelength values


elseif strcmp(file_name, 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat')==true

    % Use the TSIS-1 hypbring solar reference spectrum (mW/m^2/nm)
    % These data have 0.1 nm sampling resolution


    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [202, 2730];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of atlas_lus_modtran.txt. Must be between [200, 800] nm.', newline])
    end

    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;


elseif strcmp(file_name, 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat')==true

    % Use the TSIS-1 hypbring solar reference spectrum (mW/m^2/nm)
    % These data have 0.025 nm sampling resolution

    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [202, 2730];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of atlas_lus_modtran.txt. Must be between [200, 800] nm.', newline])
    end

    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;



elseif strcmp(file_name, 'hybrid_reference_spectrum_p025nm_resolution_c2022-11-30_with_unc.dat')==true

    % Use the TSIS-1 hypbring solar reference spectrum (mW/m^2/nm)

    % Lets check to make sure the wavelength input is within bounds of the
    % file selected

    wavelength_regime = [202, 2730];            % nanometers - wavelength boundaries

    if wavelength_boundaries(1)<wavelength_regime(1) || wavelength_boundaries(2)>wavelength_regime(2)
        error([newline, 'Wavelength is out of the range of atlas_lus_modtran.txt. Must be between [200, 800] nm.', newline])
    end

    % ------------------------------------------------------
    % -------- Reading .dat file using textscan ------------
    % ------------------------------------------------------
    % Or we could use the textscan() function instead, which allows us to define comments to ignore

    file_id = fopen([solar_source_folder,file_name], 'r');   % 'r' tells the function to open the file for reading


    format_spec = '%f %f';                                  % two floating point numbers
    B = textscan(file_id, format_spec, 'CommentStyle','#');

    index_wavelength = B{1}>=wavelength_boundaries(1) & B{1}<=wavelength_boundaries(2);

    wavelength = B{1}(index_wavelength);                % wavelengths within the user specified range

    solar_flux = B{2}(index_wavelength);                % milli-Watts/m^2/nm - flux values at the corresponding wavelength values

    % lets convert solar flux to Watts/nm/m^2

    solar_flux = solar_flux./1000;


else

    error([newline,'I dont recognize the source file you entered!',newline])

end


    % --- You must close the file! Too many open files can crash Matlab ---
    fclose(file_id);




end