%% --- Read Output Files from MIE ---
% ======================================
% The purpose of this script is to read the data output from the output
% file created by mie. The code will determine the headers for each
% column of data, and the numerical values. Headers will be stored in the
% first row of a cell array. The data will be stored in the second row of
% the cell array.


% TO-DO:
%   1) The outputs are different if we choose to output transmittance or
%   transmissivity. What are the headers?

%   2) Test to see if the new structure writing works with multiple files.



% --- By Andrew J. Buggee ---
%% --- Read in Files ---

function [dataStruct,headers,num_radii] = readMIE(path,fileName)

% How many files do we need to read?

if iscell(fileName)==true
    numFiles2Read = length(fileName);
elseif ischar(fileName)==true
    numFiles2Read = 1;
else
    error('I dont understand the input file')
end

% [fileName,path] = uigetfile('.txt');
delimeter = ' ';
headerLine = 0; % 0 if no header

data = [];

if numFiles2Read==1
    
    % old way of reading the data
%     data = importdata([path,fileName,'.OUT'],delimeter,headerLine);
    

    % new way of readng the data. cuts the time in half
    
    data_char = fileread([path,fileName,'.OUT']);
    data = cell2mat(textscan(data_char, '%f %f %f %f %f %f %f %f'));
    

elseif numFiles2Read>1
    
    for ii = 1:numFiles2Read
        
        data = cat(3,data,importdata([path,fileName{ii},'.OUT'],delimeter,headerLine));
        
    end
    
else
    
    error('There is something wrong with the number of files im supposed to read')
    
end

%% -- If there is a size distribution, we must remove the fit values --

% There is a non-homogenous size distribution if there are nan values in
% the final column of our data matrix

indexNan = isnan(data(:,end));

if sum(indexNan)>0
    
    % remove rows with nans
    data(indexNan,:) = [];
    
else
    
end


%% --- Converstion to a Structure ---


% we want an array ouput of the irradiance and radiance data for wasy
% analysis, and we want a structure output for easy manipulation on the fly



headers = cell(1,8);
headers{1} = 'wavelength';
headers{2} = 'r_eff';
headers{3} = 'refrac_real';
headers{4} = 'refrac_imag';
headers{5} = 'Qext';
headers{6} = 'ssa';
headers{7} = 'asymParam';
headers{8} = 'Qsca';





% and now we'll create the structure
% we'll start by creating the wavelength vector

% for the 'data' array, the number of rows is equal to the number of radii
% times the number of wavelengths

% Lets make sure we know how many different radii have been calculated
index = data(1,1)==data(:,1); % how many repeats of the first wavelength are there?

num_radii = sum(index); % this is the number of different radii computed by the mie code


% parse through the data

 for rr = 1:num_radii
            
       index_radii = (rr:num_radii:size(data,1))'; % index to pull out calculations with a constant radius
         
       if rr ==1
           wavelength = data(index_radii,1); % grab the wavelength values - they vary across rows
           r_eff = data(1:num_radii,2)';         % grab the values for effective radii - they varry across columns
           refrac_real = data(index_radii,3); % grab the real part of the refractive index
           refrac_imag = data(index_radii,4); % grab the imagniary part of the refractive index
       end
       
       qext(:,rr) = data(index_radii,5); % grab the extinction efficiency
       ssa(:,rr) = data(index_radii,6); % grab the single scattering albedo
       asy(:,rr) = data(index_radii,7); % grab the assymmetry parameter
       qsca(:,rr) = data(index_radii,8); % grab the assymmetry parameter
       
 end
       

if numFiles2Read==1
    % if we have to zenith view angles, then only irradiance is
    % calculated
    
    if size(data,2) == 6
        
       
        
        dataStruct = struct(headers{1},wavelength, headers{2}, r_eff,...
            headers{3},refrac_real,headers{4},refrac_imag,...
            headers{5},qext,headers{6},ssa); 
        
    elseif size(data,2) == 7
        
        dataStruct = struct(headers{1},wavelength, headers{2}, r_eff,...
            headers{3},refrac_real,headers{4},refrac_imag,...
            headers{5},qext,headers{6},ssa,...
            headers{7},asy); 
        
    elseif size(data,2) == 8
        
        dataStruct = struct(headers{1},wavelength, headers{2}, r_eff,...
            headers{3},refrac_real,headers{4},refrac_imag,...
            headers{5},qext,headers{6},ssa,...
            headers{7},asy,headers{8},qsca); 
        
        
    end
    
    
    
    
    
    
    
    
elseif numFiles2Read>1
    dataStruct = struct(headers{1,1},reshape(data(:,1,:),[],numFiles2Read),...
        headers{1,2},reshape(data(:,2,:),[],numFiles2Read),headers{1,3},...
        reshape(data(:,3,:),[],numFiles2Read),headers{1,4},...
        reshape(data(:,4,:),[],numFiles2Read),headers{1,5},...
        reshape(data(:,5,:),[],numFiles2Read),headers{1,6},...
        reshape(data(:,6,:),[],numFiles2Read),headers{1,7},...
        reshape(data(:,7,:),[],numFiles2Read),headers{1,8},...
        reshape(data(:,8,:),[],numFiles2Read));
end








end

