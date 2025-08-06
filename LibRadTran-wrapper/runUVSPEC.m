%% --- Run Intput Files with UVSPEC ---
% ======================================
% The purpose of this script is to run uvspec with the input files specfied
% by the user. The user will be required to provide to folder location of
% the input file, and the input file name. The user will also need to
% provide the output file name. This output file will be saved in the
% same folder as the input file. The input file will be fed into the
% command line in order to run uvspec.



% --- By Andrew J. Buggee ---
%% Creating the .INP file

function [inputSettings] = runUVSPEC(folderName_INP_OUT,inputName,outputName, computer_name)
%% ---- A Few Checks are Needed ----

if iscell(inputName)==true && iscell(outputName)==false
    error('inputName is a cell array, while outputName is not')

elseif iscell(inputName)==false && iscell(outputName)==true
    error('outputName is a cell array, while inputName is not')

elseif iscell(inputName)==true && iscell(outputName)==true
    if length(inputName)~=length(outputName)
        error('The number of input files doesns equal the number of output files')
    end
end

% -------------------------------------------







%% ----- Lets Read the input file -----


% Lets determine the input settings

% First we need to determine how many files we need to run
if iscell(inputName)==true
    numFiles2Run = length(inputName);
elseif ischar(inputName)==true
    numFiles2Run = 1;
else
    error('I dont understand the input file')
end



if numFiles2Run==1

    if iscell(inputName)==true
        textFile = fileread([folderName_INP_OUT,inputName{:}]);
    elseif ischar(inputName)==true
        textFile = fileread([folderName_INP_OUT,inputName]);
    end

    expr1 = '[^\n]*rte_solver [^\n]*';
    expr2 = '[^\n]*umu [^\n]*';
    expr3 = '[^\n]*phi [^\n]*';
    expr4 = '[^\n]*sza [^\n]*';
    expr5 = '[^\n]*phi0 [^\n]*';
    expr6 = '[^\n]*zout [^\n]*';
    expr7 = '[^\n]*source [^\n]*';
    expr8 = '[^\n]*wavelength [^\n]*';

    match1 = regexp(textFile,expr1,'match'); % find rte_solver typ
    match2 = regexp(textFile,expr2,'match'); % find cosine of viewing angle vector
    match3 = regexp(textFile,expr3,'match'); % find azimuth viewing angle vector
    match4 = regexp(textFile,expr4,'match'); % find the solar zenith angle
    match5 = regexp(textFile,expr5,'match'); % find the solar azimuth angle
    match6 = regexp(textFile,expr6,'match'); % find the sensor altitude
    match7 = regexp(textFile,expr7,'match'); % find the source file
    match8 = regexp(textFile,expr8,'match'); % find the wavelength range in order to trim the source file

    index1_space1 = regexp(match1{1},'\s[a-z]+'); % find the spaces
    index1_space2 = regexp(match1{1},'[a-zA-Z_0-9]\s+'); % the second index should be the last letter in the solver type

    index2_space1 = regexp(match2{1},'\s[0123456789]+'); % find a space that is followed by a number
    index2_space2 = regexp(match2{1},'[0123456789]\s+'); % find a space that comes after a number

    index3_space1 = regexp(match3{1},'\s[0123456789]+'); % find a space that is followed by a number
    index3_space2 = regexp(match3{1},'[0123456789]\s+'); % find a space that comes after a number

    index4_space1 = regexp(match4{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
    index4_space2 = regexp(match4{1},'[0123456789]\s+'); % find a space that comes after a number

    index5_space1 = regexp(match5{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
    index5_space2 = regexp(match5{1},'[0123456789]\s+'); % find a space that comes after a number - but since the varaible name, phi0, ends in a 0, we have to manually choose the second space

    index6_space1 = regexp(match6{1},'\s[0123456789]+'); % Find a space that is followed by a number
    index6_space2 = regexp(match6{1},'[0123456789]\s+'); % find a space that comes after a number
    index6_toa = regexp(match6{1}, 'toa', 'once');
    index6_allLevels = regexp(match6{1}, 'all_levels', 'once');
    index6_modelLevels = regexp(match6{1}, 'model_levels');
    index6_modelLayers = regexp(match6{1}, 'model_layers');
    index6_modelLevels_andLayers = regexp(match6{1}, 'model_levels_and_layers');

    % don't let the lack of a source stop you!
    if isempty(match7)==true

    else
        index7_space1 = regexp(match7{1},'\s[a-z]'); % find the spaces
        index7_space2 = regexp(match7{1},'[a-z]\s'); % Brackets treat the symbol literally. number of decimals tells us how many values there are in the vector
        index7_file1 = regexp(match7{1},'[/][a-z]'); % find the locaition a letter follows two dots and a forward slash
        index7_file2 = regexp(match7{1},'[.]dat');
    end

    % don't let a lack of wavelengths stop you!
    if isempty(match8)==false
        index8_space1 = regexp(match8{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index8_space2 = regexp(match8{1},'[0123456789]\s+'); % find a space that comes after a number
    end


    % determine the rte_solver type
    rte_solver = match1{1}(index1_space1(1)+1:index1_space2(2));

    % determine the umu vector
    umuStr = cell(1,length(index2_space1));

    for ii = 1:length(index2_space1)
        umuStr{ii} = match2{1}(index2_space1(ii)+1:index2_space2(ii));
    end

    umuVec = str2double(umuStr);

    % determine the phi vector
    phiStr = cell(1,length(index3_space1));

    for ii = 1:length(index3_space1)
        phiStr{ii} = match3{1}(index3_space1(ii)+1:index3_space2(ii));
    end

    phiVec = str2double(phiStr);

    % find the solar zenith angle
    sza = match4{1}(index4_space1+1:index4_space2);
    sza = str2double(sza);

    % find the solar azimuth angle
    saz = match5{1}(index5_space1+1:index5_space2(2));
    saz = str2double(saz);


    % find the sensor altitude (might be a vector!)
    zout_str = cell(1,length(index6_space1));
    zout = zeros(1,length(index6_space1));

    if isempty(zout)==false

        for nn = 1:length(zout)

            if isempty(index6_space1(nn))==false

                zout_str{nn} = match6{1}(index6_space1(nn)+1:index6_space2(nn)); % this would be for a numeric value
                zout(nn) = str2double(zout_str{nn});

            elseif isempty(index6_space1(nn))==true
                indexString1 = regexp(match6{1},'\s[a-z][a-z][a-z]'); % this would be for a string value like toa or boa
                indexString2 = regexp(match6{1},'[a-z]\s');
                zout_str{nn} = match6{1}(indexString1(1)+1:indexString2(2));

                if strcmp(zout_str,'toa')==true
                    zout(nn) = 100;
                elseif strcmp(zout_str,'sur')==true
                    zout(nn) = 0;
                end

            end

        end

    end

    if isempty(index6_toa)==false

        zout_str{end+1} = 'toa';
        zout(end+1) = 100;

    end

    if isempty(index6_allLevels)==false

        zout_str{end+1} = 'all_levels';
        zout(end+1) = inf;

    end

    if isempty(index6_modelLayers)==false

        zout_str{end+1} = 'model_layers';
        zout(end+1) = inf;

    end

    if isempty(index6_modelLevels)==false

        zout_str{end+1} = 'model_levels';
        zout(end+1) = inf;

    end

    if isempty(index6_modelLevels_andLayers)==false

        zout_str{end+1} = 'model_levels_and_layers';
        zout(end+1) = inf;

    end









    % *********************************************************
    % -------------** READING THE SOURCE FILE **---------------
    % *********************************************************

    % find the wavelength range of the output file
    if isempty(match8)==false
        wavelength_str= cell(1,length(index8_space1));

        for ii = 1:length(index8_space1)
            wavelength_str{ii} = match8{1}(index8_space1(ii)+1:index8_space2(ii));
        end
        wavelength = str2double(wavelength_str);

    else
        wavelength = [];
    end

    % ---------------------------------------------------------
    % ------ determine the solar or thermal source file -------
    % we want to store the source flux as a vector
    % all solar source files will be located in the folder: /Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/solar_flux
    % all thermal source files will be located in the foler:

    if isempty(match7)==false
        if strcmp('solar', match7{1}(index7_space1(1)+1:index7_space2(2)))

            if isempty(index7_file2)==true
                % This happens when the input is simple 'source solar' with no
                % specified file
                fileSolar = 'internal';
                source = [];
                source_flux = [];
                source_wavelength = [];

            else
                % we have defined a specfifc solar file. Figure out which
                % one it is

                if isempty(index7_file1)==false
                    % then we have a file extension with backslashes
                    fileSolar = match7{1}(index7_file1(end)+1:index7_file2(1)+3);

                else
                    % Then we have just a file name
                    fileSolar = match7{1}(index7_space1(2)+1:index7_file2(1)+3);

                end
                % Read the solar flux file over the wavelength range specified
                [source_flux, source_wavelength] = read_solar_flux_file(wavelength, fileSolar);   % W/nm/m^2
            end

        end

    else

        source = [];

    end



    % Pull all input settings into a cell array
    % first lets give them headers and labels:

    inputSettings{1,1} = 'Solver Type';
    inputSettings{1,2} = 'Cos(zva)';
    inputSettings{1,3} = 'Azimuthal Angle';
    inputSettings{1,4} = 'Solar Zenith Angle';
    inputSettings{1,5} = 'Solar Azimuthal Angle';
    inputSettings{1,6} = 'Sensor Altitude (km)';
    inputSettings{1,7} = 'Source Wavelength (nm) and Irradiance (W/nm/m^2)';

    inputSettings{2,1} = rte_solver;
    inputSettings{2,2} = umuVec;
    inputSettings{2,3} = phiVec;
    inputSettings{2,4} = sza;
    inputSettings{2,5} = saz;
    inputSettings{2,6} = zout;
    inputSettings{2,7} = [source_wavelength, source_flux];


    % if using montecarlo solver, read the mc_basename designation
    if strcmp(rte_solver, 'montecarlo')==true

        inputSettings{1,8} = 'MC Basename';

        expr_mc = '[^\n]*mc_basename [^\n]*';

        match_mc = regexp(textFile,expr_mc,'match'); % find rte_solver typ

        index_slash = regexp(match_mc{1},'[/]\w*'); % find the last forward slash

        index_lastCharacterBeforeComment = regexp(match_mc{1},'\w\s*#'); % find the last forward slash

        mc_basename = match_mc{1}(index_slash(end)+1:index_lastCharacterBeforeComment);

        inputSettings{2,8} = mc_basename;

    end



elseif numFiles2Run>1

    inputSettings = cell(numFiles2Run+1,7);

    inputSettings{1,1} = 'Solver Type';
    inputSettings{1,2} = 'Cos(zva)';
    inputSettings{1,3} = 'Azimuthal Angle';
    inputSettings{1,4} = 'Solar Zenith Angle';
    inputSettings{1,5} = 'Solar Azimuthal Angle';
    inputSettings{1,6} = 'Sensor Altitude (km)';
    inputSettings{1,7} = 'Source Wavelength (nm) and Irradiance';

    for jj=1:numFiles2Run

        textFile = fileread([folderName_INP_OUT,inputName{jj}]);

        expr1 = '[^\n]*rte_solver [^\n]*';
        expr2 = '[^\n]*umu [^\n]*';
        expr3 = '[^\n]*phi [^\n]*';
        expr4 = '[^\n]*sza [^\n]*';
        expr5 = '[^\n]*phi0 [^\n]*';
        expr6 = '[^\n]*zout [^\n]*';
        expr7 = '[^\n]*source [^\n]*';
        expr8 = '[^\n]*wavelength [^\n]*';

        match1 = regexp(textFile,expr1,'match'); % find rte_solver typ
        match2 = regexp(textFile,expr2,'match'); % find consine of viewing angle vector
        match3 = regexp(textFile,expr3,'match'); % find azimuth viewing angle vector
        match4 = regexp(textFile,expr4,'match'); % find the solar zenith angle
        match5 = regexp(textFile,expr5,'match'); % find the solar azimuth angle
        match6 = regexp(textFile,expr6,'match'); % find the sensor altitude
        match7 = regexp(textFile,expr7,'match'); % find the source file
        match8 = regexp(textFile,expr8,'match'); % find the wavelength range in order to trim the source file

        index1_space1 = regexp(match1{1},'\s[a-z]+'); % find the spaces
        index1_space2 = regexp(match1{1},'[a-z]\s+'); % the second index should be the last letter in the solver type

        index2_space1 = regexp(match2{1},'\s[0123456789]+'); % find a space that is followed by a number
        index2_space2 = regexp(match2{1},'[0123456789]\s+'); % find a space that comes after a number

        index3_space1 = regexp(match3{1},'\s[0123456789]+'); % find a space that is followed by a number
        index3_space2 = regexp(match3{1},'[0123456789]\s+'); % find a space that comes after a number

        index4_space1 = regexp(match4{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index4_space2 = regexp(match4{1},'[0123456789]\s+'); % find a space that comes after a number

        index5_space1 = regexp(match5{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index5_space2 = regexp(match5{1},'[0123456789]\s+'); % find a space that comes after a number

        index6_space1 = regexp(match6{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index6_space2 = regexp(match6{1},'[0123456789]\s+'); % find a space that comes after a number

        index7_space1 = regexp(match7{1},'\s[a-z]'); % find the spaces
        index7_space2 = regexp(match7{1},'[a-z]\s'); % Brackets treat the symbol literally. number of decimals tells us how many values there are in the vector
        %index7_file1 = regexp(match7{1},'flux[/][a-z]'); % find the location a letter follows two dots and a forward slash
        %index7_file1 = regexp(match7{1},'flux[/][a-z]'); % find the location a letter follows two dots and a forward slash
        index7_file2 = regexp(match7{1},'[.]dat');

        index8_space1 = regexp(match8{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index8_space2 = regexp(match8{1},'[0123456789]\s+'); % find a space that comes after a number


        % determine the rte_solver type
        rte_solver = match1{1}(index1_space1(1)+1:index1_space2(2));

        % determine the umu vector
        umuStr = cell(1,length(index2_space1));

        for ii = 1:length(index2_space1)
            umuStr{ii} = match2{1}(index2_space1(ii)+1:index2_space2(ii));
        end

        umuVec = str2double(umuStr);

        % determine the phi vector
        phiStr = cell(1,length(index3_space1));

        for ii = 1:length(index3_space1)
            phiStr{ii} = match3{1}(index3_space1(ii)+1:index3_space2(ii));
        end

        phiVec = str2double(phiStr);

        % find the solar zenith angle
        sza = match4{1}(index4_space1+1:index4_space2);
        sza = str2double(sza);

        % find the solar azimuth angle
        saz = match5{1}(index5_space1+1:index5_space2(2));
        saz = str2double(saz);

        % find the sensor altitude


        if isempty(index6_space1)==false
            zout = match6{1}(index6_space1(1)+1:index6_space2(2)-1); % this would be for a numeric value
            zout = str2double(zout);
        elseif isempty(index6_space1)==true
            indexString1 = regexp(match6{1},'\s[a-z][a-z][a-z]'); % this would be for a string value like toa or boa
            indexString2 = regexp(match6{1},'[a-z]\s');
            zout_str = match6{1}(indexString1(1)+1:indexString2(2));

            if strcmp(zout_str,'toa')==true
                zout = 100;
            elseif strcmp(zout_str,'boa')==true
                zout = 0;
            end
        end


        % *********************************************************
        % ---** READING THE SOURCE FILE OUTSIDE THIS FUNCTION ---**
        % *********************************************************

        % find the wavelength range of the output file
        wavelength_str = cell(1,length(index8_space1));

        for ii = 1:length(index8_space1)
            wavelength_str{ii} = match8{1}(index8_space1(ii)+1:index8_space2(ii));
        end
        wavelength = str2double(wavelength_str);


        % determine the solar or thermal source file
        % we want to store the source flux as a vector
        % all solar source files will be located in the folder: /Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/solar_flux
        % all thermal source files will be located in the foler:
        if strcmp('solar',match7{1}(index7_space1(1)+1:index7_space2(2)))

            %fileSolar = match7{1}(index7_file1(1)+5:index7_file2(1)+3);
            fileSolar = match7{1}(index7_space1(2)+1:index7_file2(1)+3);

            % Read the solar flux file over the wavelength range specified
            [source_flux, source_wavelength] = read_solar_flux_file(wavelength, fileSolar); % - (W/nm/m^2) -

        end



        % Pull all input settings into a cell array
        inputSettings{jj+1,1} = rte_solver;
        inputSettings{jj+1,2} = umuVec;
        inputSettings{jj+1,3} = phiVec;
        inputSettings{jj+1,4} = sza;
        inputSettings{jj+1,5} = saz;
        inputSettings{jj+1,6} = zout;
        inputSettings{jj+1,7} = [source_wavelength, source_flux];


        % if using montecarlo solver, read the mc_basename designation
        if strcmp(rte_solver, 'montecarlo')==true

            inputSettings{1,8} = 'MC Basename';

            expr_mc = '[^\n]*mc_basename [^\n]*';

            match_mc = regexp(textFile,expr_mc,'match'); % find rte_solver typ

            index_slash = regexp(match_mc{1},'[/]\w*'); % find the last forward slash

            index_lastCharacterBeforeComment = regexp(match_mc{1},'\w\s*#'); % find the last forward slash

            mc_basename = match_mc{1}(index_slash(end)+1:index_lastCharacterBeforeComment);

            inputSettings{jj+1,8} = mc_basename;

        end



    end
end










%% Running the .INP file from the command line

% to run a .INP file from the command line, we run uvspec by pointing the
% command line to its full location. AND the command line directory must be
% in the folder where the file you are running sits. For example, if the
% file you wish to run lives in the folder: directory/folder/runMe.file
% and the program runnning the file is in the folder:
% directory/program/fancyProgram then the code in the command line to run
% this file is:
%
% cd directory/folder/
% directory/program/fancyProgram < runMe.file > output.OUT





% --- Now we Can Run the Files ----


% --- Point to locaiton of uvspec program ---

% To run uvspec in the command line we have to point to its full location.
% To do this we will check to see what computer we are using



if strcmp('anbu8374',computer_name)

    uvspec_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/';

elseif strcmp('andrewbuggee',computer_name)

    uvspec_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/'...
        'Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/bin/'];

elseif strcmp('curc', computer_name)
    % location of the mie program
    uvspec_folderName = '/projects/anbu8374/software/libRadtran-2.0.5/bin/';

    % if running on the CURC supercomputer, you need to load the modules
    % everytime you run a job. That's because each time you run the
    % function 'system', a new unique terminal window is open. Each time a
    % new terminal window is open, the modules need to be loaded.
    cmnd_modules = ['ml purge', newline, 'ml gcc/11.2.0', newline,...
        'ml netcdf/4.8.1', newline, 'ml perl/5.36.0', newline, 'ml texlive/2021',...
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/:$PATH',...         % add libRadtran to the path
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/data/:$PATH',...    % add libRadtran data files to the path
        newline, 'export PATH=/projects/$USER/software/libRadtran-2.0.5/bin/:$PATH',...     % add uvspec location to the path
        newline, 'export PATH=', folderName_INP_OUT,':$PATH',...                            % add inp/out file locations to the path
        newline, 'export GSL_BIN=/projects/$USER/software/gsl-2.6/bin',...                  % define the binary file location for the GSL pacakge for libRadtran to find
        newline, 'export GSL_LIB=/projects/$USER/software/gsl-2.6/lib',...                  % define the library folder for the GSL package for libRadtran to fine
        newline, 'export GSL_INC=/projects/$USER/software/gsl-2.6/include',...              % define pointer for libRadtran to find the GSL pacakge
        newline, 'export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH',...                     % define the install library path for libRadtran
        newline, 'export INSTALL_DIR=/projects/$USER/software/libRadtran-2.0.5',...         % define the install directory
        newline, 'export PATH=$GSL_BIN:$PATH'];                                             % add the GSL binary files location to the path

end


% using the function 'system' runs commans in the terminal window
cmnd1 = ['cd ', uvspec_folderName];



if numFiles2Run==1

    if ischar(inputName)==true
        % cmnd2 = [uvspec_folderName,'uvspec ',...
        %            '< ',inputName,' > ', outputName];

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ',folderName_INP_OUT,inputName,' > ', folderName_INP_OUT, outputName,'.OUT',...
            ')>& ', folderName_INP_OUT,'errMsg.txt'];
        % a successful command will return a status of 0
        % an unsuccessful command will return a status of 1

    elseif iscell(inputName)==true

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ', folderName_INP_OUT, inputName{1},' > ', folderName_INP_OUT, outputName{1},'.OUT',...
            ')>&', folderName_INP_OUT, 'errMsg.txt'];

    else

        error([newline, 'I dont understand the INP filename structure', newline])

    end


    % run all commands in the terminal window
    if strcmp('curc', computer_name)

        [status] = system([cmnd_modules, ' ; ', cmnd1, ' ; ', cmnd2]);
        %[status] = system([cmnd_modules, ' ; ', cmnd2]);

    else
        [status] = system([cmnd1, ' ; ', cmnd2]);
    end

    if status ~= 0
        error(['Status returned value of ',num2str(status)])
    end



elseif numFiles2Run>1


    for ii = 1:numFiles2Run


        % cmnd2 = [uvspec_folderName,'uvspec ',...
        %            '< ',inputName,' > ', outputName];

        cmnd2 = ['(',uvspec_folderName,'uvspec ',...
            '< ', folderName_INP_OUT, inputName{ii},' > ', folderName_INP_OUT,...
            outputName{ii},'.OUT',')>& errMsg.txt'];
        % a successful command will return a status of 0
        % an unsuccessful command will return a status of 1

        if strcmp('curc', computer_name)

            [status] = system([cmnd_modules, ' ; ', cmnd1, ' ; ', cmnd2]);

        else
            [status] = system([cmnd1, ' ; ', cmnd2]);
        end


        if status ~= 0
            error(['Status returned value of ',num2str(status)])
        end
    end

else

    error('I Dont understand the file names you are trying to run')

end



end






