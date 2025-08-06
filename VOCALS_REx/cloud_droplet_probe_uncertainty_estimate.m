%% Estimate the uncertainty of the effective radius based on the measurement uncertainty of
% Droplet Measurement Technologies Cloud Droplet Probe (Lance et al.
% (2010)).


% By Andrew John Buggee

%%


function [droplet_radius_uncertainty] = cloud_droplet_probe_uncertainty_estimate(effective_radius)



% -------------------------------------------------------------------
% ---------- Uncertainty estimate based on CDP bin sizing -----------
% -------------------------------------------------------------------

% One way to define the uncertainty is using the bin spacing of the
% instrument. The bin spacing is 1 micron for bins less than 14
% microns diameter, and 2 microns for bins greater than 14 microns
% diameter. So for measured radii from 1 to 7 microns, there is an
% uncertainty of 1 micron? Well, if the bin is between 6 and 7 microns, the
% droplet could be anything between 6 and 7 microns. It's best to describe
% this bin as measuring radii of 6.5 microns, and defining the uncertainty
% as +/- 0.5 microns.


% -------------------------------------------------------------------
% ------- Uncertainty estimate based on Faber et al. (2018) ---------
% -------------------------------------------------------------------

% or we could use the results from S.Faber et al. (2018). Table 2 lists the
% absolute difference between the volume-mean diameter measured by CDP and the
% volume-mean diameter of the droplets measured with glare images of the droplets.

% VMD_difference_faber = [-1.1, 0.2, -0.1, 1.3, 1.0, 1.2, 1.5];          % microns
% mean_true_diameter = [9, 17, 24, 29, 34, 38, 46];                      % microns
%
% % The volume mean radius is equal to half the volume mean diameter
%
% VMR_difference_faber = abs(VMD_difference_faber)./2;                 % microns
% mean_true_radius = mean_true_diameter./2;                       % microns
%
%
% % if the effective radius value by vocals-rex is below 4.5 microns, assume
% % the same uncertainty as the value at 4.5 microns. If the effective radius
% % is greater than 4.5 and less than 23 microns, interpolate.
%
% droplet_radius_uncertainty = zeros(1, length(effective_radius));
%
% for nn = 1:length(effective_radius)
%
%     if effective_radius(nn)<=mean_true_radius(1)
%
%         droplet_radius_uncertainty(nn) = VMR_difference_faber(1);
%
%     elseif effective_radius(nn)>mean_true_radius(1) && effective_radius(nn)<mean_true_radius(end)
%
%         % then we interpolate
%         droplet_radius_uncertainty(nn) = interp1(mean_true_radius, VMR_difference_faber, effective_radius(nn), 'linear');
%
%     elseif effective_radius(nn)>mean_true_radius(end)
%
%         % set the uncertainty as the value of the largest radius from the
%         % Faber measurements
%         droplet_radius_uncertainty(nn) = VMR_difference_faber(end);
%
%     end
%
%
%
%
% end



% -------------------------------------------------------------------
% Uncertainty estimate based on Faber et al (2018) and Feingold (2006)
% -------------------------------------------------------------------

% CDP sizing uncertainty is quite hard to pin down. But S. Faber et al.
% (2018) state that the uncertinaty comparing higher order moments of the
% measured and true distributions are typically 10% or less. We can use the
% Feingold method of simply fudging this. Droplet diameters less than 9
% microns will have an uncertainy of 15%, and everything greater than that
% will have an uncertainty of 10%.

re_uncertainty_0_5 = 0.20;      % percentage
re_uncertainty_5_10 = 0.15;      % percentage
re_uncertainty_10_30 = 0.10;      % percentage
%re_uncertainty_all = 0.1;          % percentage


droplet_radius_uncertainty = zeros(1, length(effective_radius));

for nn = 1:length(effective_radius)

    if effective_radius(nn)<=5

        droplet_radius_uncertainty(nn) = effective_radius(nn)*re_uncertainty_0_5;

    elseif effective_radius(nn)>5 && effective_radius(nn)<=10

        % then we interpolate
        droplet_radius_uncertainty(nn) = effective_radius(nn)*re_uncertainty_5_10;

    elseif effective_radius(nn)>10

        % then we interpolate
        droplet_radius_uncertainty(nn) = effective_radius(nn)*re_uncertainty_10_30;



    end




end



% -------------------------------------------------------------------
% ------- Uncertainty estimate based on Lance et al. (2010) ---------
% -------------------- Results from figure 5 ----------------------
% -------------------------------------------------------------------

% volume mean radius uncertainty as compared to the true mean volume radius
% using laboratory experiments
% re = [4, 4.5, 5.5, 6.5, 7.5, 8.5, 19, 22];                                % microns
% re_uncertainty_boundaries = [25, 5];                   % percent uncertainty




% -------------------------------------------------------------------
% ------- Uncertainty estimate based on Lance et al. (2010) ---------
% -------------------- Results from figure 12 ----------------------
% -------------------------------------------------------------------

% radius uncertainty is a function of the droplet number concentration
% and the radius itself

% number_concentration = [10, 100, 200, 300, 400];            % #-droplets/cm^3
% re = [2.5, 3.75, 5, 10];                                % microns
% re_uncertainty = [5, 20, 40, 50, 60;...
%                   -1, 10, 20, 30, 38;...
%                   -5, 5, 10, 20, 22;...
%                   -5, -1, 5, 10, 12];                   % percent uncertainty
%
% % 2D interpolation to get the CDP uncertainty
% [Nc, Re] = meshgrid(number_concentration, re);
%
% droplet_radius_uncertainty = interp2(Nc, Re, re_uncertainty, vocalsRex.Nc, vocalsRex.re, "spline");



end

