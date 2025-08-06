%% Integrate mie calculation for a single radii over a size distribution

% INPUTS:

%   (1) r_modal - the modal radius of the droplet distribution. This is
%   NOT the ratio of the third moment of the droplet distribution to
%   the second moment. This is the most common radius observed in a large
%   random sampling of n(r), the size distribution. It is used to define
%   the droplet distributions most likely observed outcome.

%   (2) dist_var - the variance of the droplet distribution. A typical
%   value used for liquid water clouds is 7

%   (3) wavelength - (nanometers) this is the wavelength range used to calculate
%   the mie properties across the size distribution. If wavelength is a
%   single value, then a monochromatic calculation is performed. If
%   wavelength is a vector the write_mie_file will determine if it's evenly
%   spaced or not

%   (4) index_of_refraction - this is the index of refraction used in the
%       scattering calculations, effectively telling the code which substance
%       to use. Options are:
%       (a) 'water' - this is an optional string input for libRadTran
%       (b) 'ice' - this is an optional string input for libRadTran
%       (c) a + bi - this numeric option should include a real and
%       imaginary component. One mie calculation will be run for each row
%       of numeric input - thus a new mie file for each substance

%   (5) size_distribution - this is a string that tells the code which type
%   of size distribution to integrate over. The options are:
%       (a) 'gamma' - gamma droplet distribution

%   (7) computer_name - a string with the name of the computer this code on
%   which this code is running

%   (6) index - the number associated with the file that is being written.
%   This is for naming purposes so that there are unique files written, but
%   only just enough.

% OUTPUTS:
% (1) ssa_avg - single scattering albedo averaged over all drop sizes
% within the distribution defined

% (2) Qe_avg - extinction efficiency averaged over all drop sizes
% within the distribution defined

% (3) g_avg - asymmetry averaged over all drop sizes
% within the distribution defined

% (4) Qs_avg - scattering efficiency averaged over all drop sizes
% within the distribution defined




% By Andrew John Buggee
%%

function [ssa_avg, Qe_avg, g_avg, Qs_avg] = average_mie_over_size_distribution(r_eff, distribution_var, wavelength,...
    index_of_refraction, distribution_type, which_computer, index)

% ---------------------------
% ----- CHECK INPUTS --------
% ---------------------------

% The length of r_eff defines the number of unique droplet distributions.
% The length of the distribution variance must be the same length as the
% effective radius


% For the remaining entries, make sure they all have the same length as the
% length of the effective radius vector
if length(r_eff)~=length(distribution_var)
    error([newline, 'The first two inputs must be the same length',newline])
end


% % Check the wavelength input. This can only be a single value, or a vector
% % with three values
% if length(wavelength)~=1 && length(wavelength)~=3
%     error([newline, 'The wavelength input can only be a single value, implying a monochromatic calculation',...
%         newline, ', or have 3 values: the starting value, the end value, and the sampling interval',newline])
% end



% ================================================


% ------------------------------------------------
% ------------ for a gamma distribution ----------
% ------------------------------------------------

if strcmp(distribution_type, 'gamma')==true


    % For some effective radius, we have defined a droplet size distribution
    % ----- for a gamma distribution -----

    ssa_avg = zeros(1,length(r_eff));
    Qe_avg = zeros(1,length(r_eff));
    g_avg = zeros(1, length(r_eff));
    Qs_avg = zeros(1,length(r_eff));




    % for each effective radius we create a radius vector that spans the
    % full size distribution. We need to calculate the extinction
    % efficiency and the single scattering albedo for each value in our
    % radius vector


    % ---------------------------------------------------------------
    % ---------------     RUN MIE CALCULATIONS    -------------------
    % ---------------------------------------------------------------

%     % define the wavelength
%     % The wavelength input is defined as follows:
%     % [wavelength_start, wavelength_end, wavelength_step].
%     if length(wavelength)==1
%         % monochromatic calculation
%         wavelength = [wavelength, wavelength, 0];          % nanometers
%         % define the wavelength vector
%         wl_vector = wavelength(1);               % nanometers
%     elseif length(wavelength)==3
%         % broadband calculation
%         wavelength = [wavelength(1), wavelength(2), wavelength(3)];          % nanometers
%         % define the wavelength vector
%         wl_vector = wavelength(1):wavelength(3):wavelength(2);               % nanometers
%     end


    % The first entry belows describes the type of droplet distribution
    % that libRadTran uses in its calculation. The second describes the
    % distribution width. If running a mono-dispersed calculation, no
    % entry for distribution width is required.

    % For now, we run mono-dispersed calculations, and manually integrate
    % over the size distribution of choice
    % *** the second entry is the width parameter, which isn't needed for a
    % monodispersed distribution ***
    size_distribution = {'mono', []};           % droplet distribution

    % What mie code should we use to compute the scattering properties?
    mie_program = 'MIEV0';               % type of mie algorithm to run

    % Do you want a long or short error file?
    err_msg_str = 'verbose';

    
    % Lets define the radius vector the spans the full size distribution
    % By doing this outside the loop, we only have to write and run a
    % single mie calculation
    % set the total number concentration to be 1
    N0 = 1;
    
    % Use the frsit value of r_eff, but make sure r_eff is non zero
    rr = 1;
    while r_eff(rr)==0
            rr = rr+1;
    end
    [n_r, r] = gamma_size_distribution_libRadTran2(r_eff(rr), distribution_var(rr), N0);



    % Define the size of the scatterer and its scattering properties
    % Assuming a pure homogenous medium composed of a single substance.
    % The radius input is defined as [r_start, r_end, r_step].
    % where r_step is the interval between radii values (used only for
    % vectors of radii). A 0 tells the code there is no step. Finally, the
    % radii values have to be in increasing order.

    % **** The above r vector is often a different length than the vector created by libRadtran ****
    % libRadtran creates a radius vector using three inputs. Sometimes the
    % created vector is not the same length as the vector r.
    % This happens due to rounding errors. Check to make sure this doesn't
    % happen

    radius_step = mean(diff(r));

    % check the length
    while numel(r(1):radius_step:r(end))~=numel(r)

        radius_step = radius_step - 1e-15;

    end

    mie_radius = [r(1), r(end), radius_step];    % microns




    % Create a mie file
    [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
        mie_radius, wavelength, size_distribution, err_msg_str, which_computer, index);

    % run the mie file
    [~] = runMIE(mie_folder,input_filename,output_filename, which_computer);

    % Read the output of the mie file
    [ds,~,~] = readMIE(mie_folder,output_filename);

    % -----------------------------------------------------------------



    % Step through each modal radius and compute the average over the size
    % distribution
    for rr = 1:length(r_eff)


%         % ---------------- IMPORTANT ASSUMPTION -------------------
%         % the modal radius is usually less than the effective radius
%         r_modal = 0.95*r_eff(rr);
% 
%         N = dist_var(rr)^(dist_var(rr)+1)/(gamma(dist_var(rr)+1) * r_modal^(dist_var(rr)+1));  % normalization constant
% 
%         n_r = N0 * N * r.^dist_var(rr) .* exp(-dist_var(rr)*r/r_modal);                        % gamma droplet distribution


        % Compute the size distribution according to the changing effective
        % rdaius at each cloud layer
        if rr~=1
            [n_r, r] = gamma_size_distribution_libRadTran2(r_eff(rr), distribution_var(rr), N0);
        end

        % step through each wavelength. At each wavelength we integrate over
        % the vector r, a range of droplet sizes

        for ww = 1:length(wavelength)

            % Average single scattering albedo over a droplet distribution

            % ----- Old calculation below -----
            %             ssa_avg(ww,rr) = trapz(r, ds.Qsca(ww,:) .* n_r)./...
            %                 trapz(r, ds.Qext(ww,:) .* n_r);

            % ----- calculation according to Hansen and Travis (1974, pg 547-549) -----
            ssa_avg(ww,rr) = trapz(r, r.^2 .* ds.Qsca(ww,:) .* n_r)./...
                trapz(r, r.^2 .* ds.Qext(ww,:) .* n_r);




            % Compute the average asymmetry parameter over a size distribution

            % ----- Old calculation below -----
            %             g_avg(ww,rr) = trapz(r, ds.asymParam(ww,:) .* n_r)./...
            %                 trapz(r, n_r);

            % ----- calculation according to Kokhanovky (Cloud optics, pg 69) -----
            g_avg(ww,rr) = trapz(r, ds.asymParam(ww,:) .* ds.Qsca(ww,:) .* r.^2 .* n_r)./...
                trapz(r, ds.Qsca(ww,:) .* r.^2 .* n_r);




            % Compute the average extinction efficiency over a droplet size
            % distribution

            % ----- Old calculation below -----
            %             Qe_avg(ww,rr) = trapz(r, ds.Qext(ww,:) .* n_r)./...
            %                 trapz(r, n_r);

            % ----- calculation according to Hansen and Travis (1974, pg 547-549) -----
            Qe_avg(ww,rr) = trapz(r, r.^2 .* ds.Qext(ww,:) .* n_r)./...
                trapz(r, r.^2 .* n_r);



            
            % Compute the average scattering efficiency over a droplet size
            % distribution


            % ----- calculation according to Hansen and Travis (1974, pg 547-549) -----
            Qs_avg(ww,rr) = trapz(r, r.^2 .* ds.Qsca(ww,:) .* n_r)./...
                trapz(r, r.^2 .* n_r);



        end

    end





else

    error([newline, 'I can only integrate over gamma distributions at the moment...',newline])

end






end