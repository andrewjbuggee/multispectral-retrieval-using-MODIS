%% Create a txt file where the second column is the non-uniformly spaced wavelength grid


% By Andrew John Buggee

%%

function wl_filename = write_nonuniformly_spaced_wavelength_file(wavelength, index)

    % write a txt file and save it to the mie_calcualtions folder

    % Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

elseif strcmp(computer_name,'andrewbuggee')==true

    mie_folder = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

elseif strcmp(computer_name,'curc')==true

    mie_folder = '/projects/anbu8374/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';


end

%%

% create the filename
wl_filename = [mie_folder, 'wavelength_grid_',num2str(index),'.txt'];

% Create the wavelength file
fileID = fopen(wl_filename, 'w');

% first column is the number associated with the index of the wavelength
% the second column is the wavelength

% make sure wavelength is a column vector
if size(wavelength,1)==1 && size(wavelength,2)>1
    wavelength = wavelength';
end

% create the wavelength indexs
wavelength_index = (1:length(wavelength))';

fprintf(fileID, '%f %f \n', wavelength_index, wavelength);


% Close the file!
fclose(fileID);









end