%% Read and Interpolate Mie Calculations

% This function will read and interpolate precompute Mie calculations for
% water droplets of varrying radii.

% INPUTS:
%   (1) xq - querry points. If one wants to know the mie properties for a
%   droplet of radius r at wavelength w, xq is a vector that defines the
%   locations r and w for interpolation. xq = [w,r] so please provide the
%   wavelength first and then the radius. The wavelength should be in
%   nanometers and the radius should be in microns. Sorry! 
%
%   (2) distribution - choose between two lookup tables:
%       (a) 'gamma' - one using a gamma distribution of droplets with an alpha
%       factor of 7
%       (b) 'mono' - one using a single particle of precisely the diameter in
%       question
%
%   (3) justQ_flag - tells the code you only need Qext - The entire
%   Mie_Properties file is 29 MB, which takes a while to load! To write 10
%   wc files, reading in the Mie_Properties file 10 times, it takes 8
%   seconds. But most of the time we only need Qext, so there is a file for
%   both the monodispersed and the gamma distribution that contains only
%   the computed values for Qext, which is just 2.9 MB.

% OUTPUTS:
%   (1) yq - the values of each mie parameter interpolated at the locations
%   specified by xq. The variables in the output matrix are:
%       (a) wavelength (nm)
%       (b) effective radius (microns)
%       (c) refrac_real - real part of the refractive index
%       (d) refrac_imag - imaginary part of the refractive index
%       (e) q_ext - the extinction efficiency
%       (f) omega - the single scattering albedo
%       (g) gg - asymmetry parameter
%       (h) q_scat - scattering efficiency

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

% ------------------------------------------------------------------------
%%

function [yq] = interp_mie_computed_tables(xq,distribution, justQ_flag)

% ----------------------------------------
% ------------- Check Inputs! ------------
% ----------------------------------------

% Check to make sure there are two inputs


if nargin~=3
    error([newline,'Not enough inputs. Need 3: the first defines the query points, ',...
        'the next chooses the droplet size distribution, and the last tells the code to load ',...
        'the entire Mie_Properties table, or just the Q_ext values.', newline])
end

% Check to make sure xq is a vector with two values

if size(xq,2)~=2
    error([newline,'The query vector xq must have two inputs, one for wavelength and one for radii', newline])
end

% Determine which computer you're using
computer_name = whatComputer;

% use the proper path
if strcmp(computer_name,'anbu8374')==true
    folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    
elseif strcmp(computer_name,'andrewbuggee')==true
    folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
    
end

%%
% ------------------------------------------
% ------------- Function Stuff! ------------
% ------------------------------------------




% The data fields are in the following order:
%       (1) wavelength (nanometers)
%       (2) effective radius (microns)
%       (3) refrac_real - real part of the refractive index
%       (4) refrac_imag - imaginary part of the refractive index
%       (5) q_ext - the extinction efficiency
%       (6) omega - the single scattering albedo
%       (7) gg - asymmetry parameter
%       (8) q_scat - scattering efficiency

num_calcs = size(xq,1);                             % This is the number of interpolation points

if strcmp(distribution, 'gamma')==true
    
    % Let's load the mie compute look-up table for a gamma droplet
    % distribution
    
    if justQ_flag==true
        
        % if this is true, we only load the Q_ext calculations!
        filename = 'Q_ext_1nm_sampling_gamma_7.txt';
        format_spec = '%f';        % 1 column of data
        
        
        
        % ----- READ IN DATA USING TEXTSCAN() ---------------
        
        file_id = fopen([folder_path,filename]);
        
        data = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
            {'\r', ' '}, 'MultipleDelimsAsOne',1);
        
        data = reshape(data{1},100,[]);                                     % rehsape the data into a matrix
        
        % -- Define the boundaries of interpolation for wavelength and radius --
        
        wavelength_bounds = [100, 3000];            % nanometers - wavelength boundaries
        r_eff_bounds = [1, 100];                    % microns - effective radius boundaries
        
        % ---- Create the grids for interpolation ----
        
        wl = wavelength_bounds(1):wavelength_bounds(2);         % nm
        r_eff = r_eff_bounds(1):r_eff_bounds(2);                % microns
        
        [WL, R_eff] = meshgrid(wl,r_eff);                       % Meshgrid for inerpolation
        
        if any(xq(:,1) < wavelength_bounds(1)) || any(xq(:,1) > wavelength_bounds(2)) || any(xq(:,2) < r_eff_bounds(1)) || any(xq(:,2) > r_eff_bounds(2))
            
            % if any of these are true, then we will extrapolate
            error(['Query points are outside the bounds of the data set. The acceptable ranges are:',newline,...
                'Wavelength: [100, 3000] nm', newline,...
                'Effective Radius: [1, 100] microns', newline]);
            
        else
            
            % then we will interpoalte
            % Lets grab all of the values we need in the data set
            
            yq = interp2(WL, R_eff, data, xq(:,1), xq(:,2));
            
            
            % Lets include the wavelength and effective radius in the
            % interpolated data cube
            yq = [xq, yq];
            
            
            
        end
        
    else
        
        filename = 'Mie_calcs_1nm_sampling_gamma_7.OUT';          % This file computes re=1:1:100 and lambda=100:1:3000
        
        %filename = 'Mie_calcs_gamma7_10nm_sampling_txt.OUT';                % This file computes re=1:10:100 and lambda=100:10:2500
        
        format_spec = '%f %f %f %f %f %f %f %f';        % 8 columns of data
        
        file_id = fopen([folder_path,filename]);
        
        data_raw = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
            {'\r', ' '}, 'MultipleDelimsAsOne',1);
        
        % Set up the zero array
        yq = zeros(num_calcs,size(data_raw,1)-2);                                           % The first two rows are not needed
        
        
        % -- Define the boundaries of interpolation for wavelength and radius --
        
        wavelength_bounds = [data_raw{1}(1), data_raw{1}(end)];            % nanometers - wavelength boundaries
        r_eff_bounds = [data_raw{2}(1), data_raw{2}(end)];                    % microns - effective radius boundaries
        
        
        % ---- Create the grids for interpolation ----
        
        wl = unique(data_raw{1});         % nm
        r_eff = unique(data_raw{2});                % microns
        
        [WL, R_eff] = meshgrid(wl,r_eff);                       % Meshgrid for inerpolation
        
        
        % ------------------------------------------------------------
        % Let's check to make sure the wavelength and effective radius
        % query locations are within the bounds of the table. If either
        % are outside the bounds of interpolation, we will extrapolate
        % and issue a warning
        % ------------------------------------------------------------
        
        if any(xq(:,1) < wavelength_bounds(1)) || any(xq(:,1) > wavelength_bounds(2)) || any(xq(:,2) < r_eff_bounds(1)) || any(xq(:,2) > r_eff_bounds(2))
            
            % if any of these are true, then we will extrapolate
            error(['Query points are outside the bounds of the data set. The acceptable ranges are:',newline,...
                'Wavelength: [100, 3000] nm', newline,...
                'Effective Radius: [1, 100] microns', newline]);
            
        else
            
            % then we will interpoalte
            % Lets grab all of the values we need in the data set
            
            for ii = 3:length(data_raw)                         % The first two values are wavelength and effective radius
                
                data = reshape(data_raw{ii}, length(unique(data_raw{2})),[]);
                yq(:,ii-2) = interp2(WL, R_eff, data, xq(:,1), xq(:,2));
                
            end
            
            
            % Lets include the wavelength and effective radius in the
            % interpolated data cube
            yq = [xq, yq];
            
            
            
        end
        
        
    end
    
    
elseif strcmp(distribution, 'mono')==true
    
    % Let's load the mie compute look-up table for a monodispersed droplet
    % distribution
    
    
    
    if justQ_flag==true
        
        % if this is true, we only load the Q_ext calculations!
        filename = 'Q_ext_1nm_sampling_monodispersed_radii_1_to_500_microns_reduced.txt';
        format_spec = '%f';        % 1 column of data
        
        
        
        % ----- READ IN DATA USING TEXTSCAN() ---------------
        
        file_id = fopen([folder_path,filename]);
        
        data = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
            {'\r', ' '}, 'MultipleDelimsAsOne',1);
        data = reshape(data{1},340,[]);                                     % rehsape the data into a matrix
        
        
        % -- Define the boundaries of interpolation for wavelength and radius --
        
        wavelength_bounds = [100, 3000];            % nanometers - wavelength boundaries
        r_eff_bounds = [1, 500];                    % microns - effective radius boundaries
        
        
        % ---- Create the grids for interpolation ----
        
        wl = wavelength_bounds(1):1:wavelength_bounds(2);         % nm
        % From the defined radii values listed in the file above
        r_eff = [r_eff_bounds(1):300, 305:5:r_eff_bounds(2)];                % microns
        
        [WL, R_eff] = meshgrid(wl,r_eff);                       % Meshgrid for inerpolation
        
        
        
        if any(xq(:,1) < wavelength_bounds(1)) || any(xq(:,1) > wavelength_bounds(2)) || any(xq(:,2) < r_eff_bounds(1)) || any(xq(:,2) > r_eff_bounds(2))
            
            % if any of these are true, then we will extrapolate
            error(['Query points are outside the bounds of the data set. The acceptable ranges are:',newline,...
                'Wavelength: [100, 3000] nm', newline,...
                'Effective Radius: [1, 100] microns', newline]);
            
        else
            
            % then we will interpoalte
            % Lets grab all of the values we need in the data set
            yq = interp2(WL, R_eff, data, xq(:,1), xq(:,2));
            
            
            % Lets include the wavelength and effective radius in the
            % interpolated data cube
            yq = [xq, yq];
            
            
            
        end
        
    else
        
        % The data file will be parsed into 8 cells of data
        %   (1) wavelength
        %   (2) effective radius
        %   (3) index of refraction real
        %   (4) index of refraction imaginary
        %   (5) Extinction efficiency
        %   (6) single scattering albedo
        %   (7) asymmetry parameter
        %   (8) scattering efficiency
        
        filename = 'Mie_calcs_monodispersed_Bohren_Huffmam_radii_1_to_500_microns_reduced_file.txt';
        format_spec = '%f %f %f %f %f %f %f %f';        % 8 columns of data
        

        file_id = fopen([folder_path,filename]);
        
        data_raw = textscan(file_id, format_spec, 'CommentStyle','#', 'Delimiter',...
            {'\r', ' '}, 'MultipleDelimsAsOne',1);
        
        % Set up the zero array
        yq = zeros(num_calcs,length(data_raw)-2);                                           % The first two rows are not needed
        
        
        % -- Define the boundaries of interpolation for wavelength and radius --
        
        wavelength_bounds = [data_raw{1}(1), data_raw{1}(end)];            % nanometers - wavelength boundaries
        r_eff_bounds = [data_raw{2}(1), data_raw{2}(end)];                    % microns - effective radius boundaries
        
        
        % ---- Create the grids for interpolation ----
        
        wl = unique(data_raw{1});         % nm
        r_eff = unique(data_raw{2});                % microns
        
        [WL, R_eff] = meshgrid(wl,r_eff);                       % Meshgrid for inerpolation
        
        
        % ------------------------------------------------------------
        % Let's check to make sure the wavelength and effective radius
        % query locations are within the bounds of the table. If either
        % are outside the bounds of interpolation, we will extrapolate
        % and issue a warning
        % ------------------------------------------------------------
        
        if any(xq(:,1)< wavelength_bounds(1)) || any(xq(:,1)>wavelength_bounds(2)) || any(xq(:,2) < r_eff_bounds(1)) || any(xq(:,2) > r_eff_bounds(2))
            
            % if any of these are true, then we will extrapolate
            error(['Query points are outside the bounds of the data set. The acceptable ranges are:',newline,...
                'Wavelength: [100, 3000] nm', newline,...
                'Effective Radius: [1, 100] microns', newline]);
            
        else
            
            % then we will interpoalte
            % Lets grab all of the values we need in the data set
            
            for ii = 3:length(data_raw)                         % The first two columns of data are wavelength and effective radius
                
                data = reshape(data_raw{ii}, length(unique(data_raw{2})),[]);               % variables vaires with effective radius along rows and wavelength along columns
                yq(:,ii-2) = interp2(WL, R_eff, data, xq(:,1), xq(:,2));
                
            end
            
            
            % Lets include the wavelength and effective radius in the
            % interpolated data cube
            yq = [xq, yq];
            
            
            
        end
        
        
    end
    
else
    error([newline,'The distribution you provided is not valid. Enter "gamma" or "mono"',newline])
    
    
    
end



end



% ----- READ IN DATA USING IMPORTDATA() ---------------


%     delimeter = ' ';
%     headerLine = 0; % 0 if no header
%     data = importdata([folder_path,filename],delimeter,headerLine);
%     data = reshape(data',8, 100, []);                                       % converts to data cube (row,col,depth) = (parameters,r_eff, lambda)
%
%
%
%
%     % Set up the zero array
%     yq = zeros(1,size(data,1)-2);                                           % The first two rows are not needed
%
%     if any(xq(:,1)<= wavelength_bounds(1)) || any(xq(:,1)>=wavelength_bounds(2)) || any(xq(:,2) <= r_eff_bounds(1)) || any(xq(:,2) >= r_eff_bounds(2))
%
%         % if any of these are true, then we will extrapolate
%         error(['Query points are outside the bounds of the data set. The acceptable ranges are:',newline,...
%             'Wavelength: [100, 3000] nm', newline,...
%             'Effective Radius: [1, 100] microns', newline]);
%
%     else
%
%         % then we will interpoalte
%         % Lets grab all of the values we need in the data set
%         for nn = 1:num_calcs
%
%             for ii = 1:(size(data,1)-2)                         % The first two values are wavelength and effective radius
%
%                 data2interpolate = reshape(data(ii+2,:,:), 100,[]);
%                 yq(nn,ii) = interp2(WL, R_eff, data2interpolate, xq(nn,1), xq(nn,2));
%
%             end
%
%
%         end
%
%         % Lets include the wavelength and effective radius in the
%         % interpolated data cube
%         yq = [xq, yq];
%
%
%
%     end


