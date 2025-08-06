%% This function will write a .INP file for a Mie calculation using libRadTran

% INPUTS:
%   (1) mie_program - the type of mie program to use to calculate
%   scattering properties. Entries can be either:
%       (a) BH - the bohren huffman algorithm
%       (b) MIEV0 - the wiscombe algorithm

%   (2) index_refraction - this is the index of refraction used in the
%   scattering calculations, effectively telling the code which substance
%   to use. Options are:
%       (a) 'water' - this is an optional string input for libRadTran
%       (b) 'ice' - this is an optional string input for libRadTran
%       (c) a + bi - this numeric option should include a real and
%       imaginary component. One mie calculation will be run for each row
%       of numeric input - thus a new mie file for each substance

%   (3) re - (units: microns) - this entry sets the radius of the sphere
%   for mie calculations. If a distribution is used (gamma or lognormal) it
%   is the modal radius of the distribution. The input is of the following
%   form: [r_start, r_end, r_step]. Where r_step defines the interval
%   spacing and therefore the number of radii to compute. If the third
%   entry is a zero, that means there is no spacing. This works for either
%   a single radii or two radii.

%   (4) wavelength - (nanometers) - wavelength of light to use in the
%   calculation. The total number of calculations per file are (number of
%   radii * number of wavelengths). This input is a vector of wavelengths,
%   with equal grid spacing or not. LibRadTran requires the input to be in 
%   the same format as the effective radius:[wavelength_start, 
%   wavelength_end, wavelength_step]. This is only true if the wavelength
%   spacing is constant. If not, you can specify a wavelength file.
%   Wavelength step happens to be a separate input in libRadTran, and is
%   not required if you're running a calculation at a single wavelength,
%   but a starting and ending wavelength is always required.

%   (5) distribution - This entry describes the droplet distribution to use
%   (first entry) and the distribution width (second entry). Types of
%   distribution entries that are accepted:
%       (a) 'mono' - monodispersed distribution
%       (b) 'gamma' - gamma droplet distribution. By default this will use
%       a gamma distribution with an alpha value of 7, which is typical for
%       liquid water clouds.
%       (c) 'lognormal' - a log normal distribution. The second input is
%       the sigma value which describes the spread


%   (6) err_msg_str - this describes how long the error message should be.
%   There are two possible inputs:
%       (a) 'quiet' - this only prints the essential information
%       (b) 'verbose' - this prints a very long and detailed message about
%       the calculation.

%   (7) computer_name - a string with the name of the computer this code on
%   which this code is running

%   (8) index - a value to append at the end of the file name so that no
%   two filenames are the same while parallel computing is going on


% By Andrew John Buggee

%%

function [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_refraction,...
            re, wavelength, distribution, err_msg_str, computer_name, index1, indexRadius)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 8 inputs


if nargin~=9
    error([newline,'Not enough inputs. Need 8: mie program type, index of refraction, droplet effective radius,',...
        [' wavelength, droplet distribution width, the error message command, the computer name,',...
         'the file number, and the radius for the file naming system'], newline])
end


% The first cell entry in "distribution" is the distribution type. The
% second cell entry is the distribution width.

if strcmp(distribution{1}, 'mono')==false && strcmp(distribution{1}, 'gamma')==false

    error([newline,'I dont recognize the droplet distribution. Must be either "mono" or "gamma"', newline])
end


if strcmp(mie_program, 'MIEV0')==false && strcmp(mie_program, 'BH')==false

    error([newline,'I dont recognize the mie program. Must be either "MIEV0" or "BH"', newline])
end




% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

elseif strcmp(computer_name,'andrewbuggee')==true

    mie_folder = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

elseif strcmp(computer_name,'curc')==true

    
    mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/';

    
    if ~exist(mie_folder, 'dir')

        mkdir(mie_folder)
    end



end

%% Determine how many files need to be written

% We can run a single mie file at different radii and different
% wavelength, but only a single index of refraction can be used at a time

if ischar(index_refraction)
    % There there is only one type of scatterer we are dealing with
    numFiles = 1;

elseif iscell(index_refraction)
    % This means there are multiple scatterers using the libRadTran string
    % entries ('water' or 'ice')
    numFiles = length(index_refraction);

elseif isnumeric(index_refraction)
    % This means the index of refraction was entered manually by the user.
    % Check to see how many indices were entered
    numFiles = length(index_refraction);

end




%%
% ------------------------------------------------------------
% ---------------------- WRITE MIE FILE -----------------------
% ------------------------------------------------------------

% Create comments for each line
comments = {'# Mie code to use', '# refractive index to use', '# specify effective radius grid (microns)',...
    '# Specify size distribution and distribution width','# Define wavelength boundaries (nanometers)',...
    '# Define interval of wavelength sampling', '# define output variables','# error file length'};


for ff = 1:numFiles

    % --------------------------------------------
    % *** create the input and output filename ***
    % --------------------------------------------

    if numFiles==1

        if ischar(index_refraction)==1

            input_filename = ['Mie_calc_refrac_',index_refraction,'_distribution_',distribution{1},...
                '_nn-',num2str(index1),'_', num2str(indexRadius),'.INP'];
            output_filename = ['OUTPUT_',input_filename(1:end-4)];

            % Create the water cloud file
            fileID = fopen([mie_folder,input_filename], 'w');

        elseif isnumeric(index_refraction)==1
            input_filename = ['Mie_calc_refrac_',num2str(index_refraction),'_distribution_',distribution{1},...
                 '_nn-',num2str(index1),'_', num2str(indexRadius),'.INP'];
            output_filename = ['OUTPUT_',input_filename(1:end-4)];

            % Create the water cloud file
            fileID = fopen([mie_folder,input_filename], 'w');

        end

    elseif numFiles>1

        if iscell(index_refraction)==1

            input_filename{ff} = ['Mie_calc_refrac_',index_refraction{ff},'_distribution_',distribution{1},...
                 '_nn-',num2str(index1),'_', num2str(indexRadius),'.INP'];
            output_filename{ff} = ['OUTPUT_',input_filename(1:end-4)];

            % Create the water cloud file
            fileID = fopen([mie_folder,input_filename{ff}], 'w');

        elseif isnumeric(index_refraction)==1
            input_filename{ff} = ['Mie_calc_refrac_',num2str(index_refraction(ff,:)),'_distribution_',...
                distribution{1}, '_nn-',num2str(index1),'_', num2str(indexRadius),'.INP'];
            output_filename{ff} = ['OUTPUT_',input_filename(1:end-4)];

            % Create the water cloud file
            fileID = fopen([mie_folder,input_filename{ff}], 'w');

        end

    end




    % fprintf writes lines in our text file from top to botom
    % .INP files for mie calculations always require the same inputs

    % to write column vectors in a text file, we have to store them as row
    % vectors

    % ----------------------------------
    % Define the mie program code to use
    % ----------------------------------

    fprintf(fileID, '%11s %5s          %s \n','mie_program',mie_program,comments{1});





    % ----------------------------------
    % Define the Index of Refraction!
    % ----------------------------------

    % check to see if the index of refraction is a string or a number
    if isstring(index_refraction)==true || ischar(index_refraction)==true
        fprintf(fileID, '%6s %s          %s \n','refrac', index_refraction, comments{2});
    elseif isnumeric(index_refraction)==true
        % if true we to tell the code we have a user input value for the index
        % of refraction. NOTE: both the real and compelx parts have to be
        % entered as positive numbers
        fprintf(fileID, '%11s %f %f          %s \n','refrac user', real(index_refraction), imag(index_refraction), comments{2});

    else

        error([newline,'I dont recognize the index of refraction input.',newline])
    end



    % ---------------------------------------------------------------------
    % Write in the value for the modal radius. Check to see if its a vector
    % ---------------------------------------------------------------------

    if length(re)==3
        % If this is true, we create an equally spaced radial vector

        % The input vector re is ordered as follows: re = [r_start, r_end,
        % r_step]

        re_start = re(1);
        re_end = re(2);
        re_step = re(3);

        % The step value has to be written with many significant digits to
        % ensure rounding doesn't prevent the full vector from being computed.
        if re_start~=re_end && re_step~=0
            fprintf(fileID,'%5s %3.8f %3.8f %3.15f          %s \n', 'r_eff', re_start, re_end, re_step, comments{3});

        elseif re_start==re_end && re_step==0

            % there is only a single value for re
            fprintf(fileID,'%5s %3.2f          %s \n', 'r_eff', re(1), comments{3});

        else

            error([newline, 'The input r_eff is an empty array', newline])

        end

    elseif length(re)==1

        % If this is true, then there is only a single radius value to compute
        fprintf(fileID,'%5s %3.5f          %s \n', 'r_eff', re, comments{3});

    else

        error([newline,'There can only be a single entry for the raidus, or 3 entries that define a vector.',newline])

    end






    % ----------------------------------------------------------------
    % Write in the value for the droplet distribution, if its not mono
    % ----------------------------------------------------------------

    if strcmp(distribution{1},'gamma')==true

        fprintf(fileID,'%12s %5s %f         %s \n', 'distribution', distribution{1}, distribution{2}, comments{4});

    elseif strcmp(distribution{1},'lognormal')==true

        fprintf(fileID,'%12s %5s %f         %s \n', 'distribution', distribution{1}, distribution{2}, comments{4});

    elseif strcmp(distribution{1},'mono')==true

        % if this is true, we don't need to write a distribution line

    else

        error([newline, 'I dont recognize the distribution string.', newline])

    end




    % ---------------------------------------
    % ----- define the wavelength range -----
    % ---------------------------------------
    
    % if wavelength has only a single entry, then this is a monochromatic
    % calcualtion
    if length(wavelength)==1

        wavelength_start = wavelength(1);           % nm
        wavelength_end = wavelength(1);             % nm

        fprintf(fileID,'%10s  %5.2f %5.2f          %s \n', 'wavelength', wavelength_start, wavelength_end, comments{5});

    % determine if wavelength is a vector and if the vector is evenly spaced    
    elseif length(wavelength)==3 && isuniform(wavelength)==false
        
        % if true, the first entry is the starting wavelength, the second
        % entry is the final wavelength, and the 3rd entry is the sampling
        % rate.

        % define the starting and ending wavelength
        wavelength_start = wavelength(1);           % nm
        wavelength_end = wavelength(2);             % nm
        % compute the grid spacing
        wavelength_step = wavelength(3);        % nm

        fprintf(fileID,'%10s  %5.2f %5.2f          %s \n', 'wavelength', wavelength_start, wavelength_end, comments{5});
        fprintf(fileID,'%15s  %5.2f          %s \n', 'wavelength_step', wavelength_step, comments{6});

    elseif length(wavelength)>1 && isuniform(wavelength)==false

        % If the wavelength vector doesn't have uniform spacing, write a
        % .txt file where the second column is the wavelength vector
        wl_filename = write_nonuniformly_spaced_wavelength_file(wavelength, index1);


        fprintf(fileID,'%10s  %s          %s \n', 'wavelength', wl_filename, comments{5});

    else

        error([newline, 'I think the wavelength vector is empty.', newline])

    end






    % ---------------------------
    % Define the output variables
    % ---------------------------

    % But first make a comment
    fprintf(fileID,'\n%s \n', comments{7});
    fprintf(fileID,'%11s  %s %s %s %s %s %s %s %s \n', 'output_user','lambda', 'r_eff', 'refrac_real',...
        'refrac_imag', 'qext', 'omega', 'gg', 'qsca');





    % -----------------------
    % Print the error message
    % -----------------------
    fprintf(fileID, '%s          %s', err_msg_str, comments{8});


    % Close the file!
    fclose(fileID);



end


end



