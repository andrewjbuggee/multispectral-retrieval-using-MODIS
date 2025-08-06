%% ----- CALCULATE THE REFLECTANCE FUNCTION FOR MODIS-LIKE Retrievals-----

% *** IMPORTANT ***  MODIS defines reflectance differently than the ATBD.
% MODIS does not define reflectance with a cosine of the solar zenith
% angle in the denominator, which is how the MODIS cloud retrieval ATBD
% defines reflectance. To properly compare reflectances with MODIS, do 
% NOT include this term in the denominator



% INPUTS:
%   (1) inputSettings

%   (2) dataStructure  - output results from uvspec

%   (3) spec_response - spectral response function


% OUTPUTS:
%   (1) R - This is the reflectance over a finite bandwidth. R integrates
%   the monochromatic reflectance, R_lambda, over the wavelength range

%   (2) R_lambda - This is the reflectance for a monochromatic calculation.
%   If the equation of radiative transfer is solved for a single
%   wavelength, then we don't need to integrate. 


% By Andrew J. Buggee
%% ------ Read input settings and output data from uv_spec -----

function [R,R_lambda] = reflectanceFunction_4modis(inputSettings,ds, spec_response)

% Geometry values from input Settings -
mu = inputSettings{2}; % cosine of the viewing zenith angle
phi = inputSettings{3}; % sensor aziumuth angle
phi0 = inputSettings{5}; % solar azimuth angle
sensorAlt = inputSettings{6}; % sensor altitude in km
source = inputSettings{7}; % - W/(m^2 nm) - source irradiance

% radiative transfer solutions
wavelength = ds.wavelength;     % nm
irrad0 = source(:,2).*1e3; % - mW/(m^2 nm) -  source irradiance


% a few other constants calculated from the inputs
geomSets = length(mu)*length(phi);



%% ----- CALCULATE REFLECTANCE FUNCTION -----

% reflectance varies with lambda, tau_cloud, r_e, viewing angle, solar
% zenith angle, and the viewing azimuth angle. uvspec can only vary the
% viewing azimuth and viewing zenith angle for a given file. So we will
% create a 3D array that varies these three parameters per file

% ***** Radiance has units of mW/m^2/nm/sr *****

R_lambda = zeros(length(wavelength),geomSets); % phi (azimuth) changes first, then mu (cos(sza))
R = zeros(length(mu),length(phi));
% if there is a single monochromatic wavelength then we don't need to
% integrate. We simply divide
if length(wavelength)==1
    
    for ii = 1:geomSets
        
        R_lambda(ii) = pi*ds.radiance(ii).value./(irrad0); % - 1/sr - reflectance function for monochromatic calculation
        
    end
    
elseif length(wavelength)>1
    
    for ii = 1:geomSets
        
        % First calculate the reflectance function at each discrete
        % wavelength within the wavelength band
        R_lambda(:,ii) = pi*ds.radiance(ii).value./(irrad0); % - 1/sr/nm - reflectance function for monochromatic calculation
        R(ii) = trapz(wavelength,R_lambda(:,ii).*spec_response.*irrad0)./trapz(wavelength,spec_response.*irrad0); % - 1/sr - reflectance function over a finite bandwidth
        
    end
    
else
    
    error('Something is wrong with the wavelength vector');
    
end









end