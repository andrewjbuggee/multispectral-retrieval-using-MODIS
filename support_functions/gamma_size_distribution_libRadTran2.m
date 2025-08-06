% Create a gamma particle size distribution using the libRadTran definition

% INPUTS:
%   (1) r_eff - effective radius (microns) - this is the ratio of the third
%   moment of the distribution to the second moment. It is used to define
%   the constant, b, which resides within the exponential

%   (2) alpha - This is the second argument in the definition of the gamma
%   distribution within libRadTran. Larger values increase the width of the
%   distribution

%   (3) N0 - total droplet concentration (# molecules / unit of volume) - 
%   this is the total droplet number concentration including all sizes.
%   If n(r) is integrated over r, we get N0.

% OUTPUTS:
%   (1) n(r) - 1/(microns * unit volue) - the number concentration of droplets 
%   for a given radius r - this output is a vector of the same length as
%   r, which is a hard-coded vector. Whatever the units of the total number
%   concentration are will be passed on to the output

%   (2) r - the independent variable to defines our distribution n(r). 


% By Andrew John Buggee

%%

function [n_r,r] = gamma_size_distribution_libRadTran2(r_eff, alpha, N0)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=3
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end



if N0<0
    
    error([newline,'The number concentration must be greater than 0!', newline])
end

if r_eff<0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end


%% Define the gamma distribution
% This definition follows the definition of the libRadTran manual, equation
% 6.10

% For some effective radius, we have define a droplet size distribution over a
% range of radii values

r = linspace(0.001*r_eff, 7*r_eff, 300);                  % microns - vector based on C.Emde (2016)


% --- first, solve for the constant, b ---
b = gamma(4 + alpha) / (r_eff * gamma(3 + alpha));

% --- Using Mathematica to solve for the normalization consant ---
N = 1/(b^(-1 - alpha) * gamma(1 + alpha));

% --- define the full distribution ---
n_r = N0 * N * r.^alpha .* exp(-b .* r);                            % 1/microns/unit-volume - gamma droplet distribution




end