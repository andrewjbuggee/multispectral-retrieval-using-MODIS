%% This function will create a cloud droplet profile based on one of four
% physical assumptions. These are taken from S.Platnick's 2000 paper

% INPUTS:
%   (1) re_topBottom - effective droplet radius (microns) - these are the
%   boundary values for the profile you wish to create. You enter them as
%   a vector in the order specified by the variable name [re_top,
%   re_bottom], or, if retrieving r_middle, enter this as [re_top,
%   re_middle, re_bottom]

%   (2) zT - the vertical independent variable, which is
%   either geometric altitude or optical depth. If using geometric altitude
%   (z), this variable should be defined in units of kilometers and should
%   start with z_bottom and progress towards z_top. If using optical depth
%   (tau), this vector should start with 0 (cloud top) and end with the
%   total cloud optical thickness.

%   (3) independent varaible - a string that tells the code what the
%   vertical indepdendent variable is. There are two possible inputs here:
%       (a) 'altitude' - tells the code to compute r as a function of
%       geometric height
%       (b) 'optical_depth' - tells the code to compute r as a function of
%       optical thickness within the cloud

%   (4) constraint - the physical constraint (string) - there are four
%   different string options for a physical constraint:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.

% OUTPUTS:
%   (1) re - effective droplet radius profile (microns). The starts at
%   cloud top and descends to the cloud bottom value


% By Andrew John Buggee
%%

function re = create_droplet_profile2(re_topBottom, zT, independentVariable, constraint)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 4 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=4
    error([newline,'Not enough inputs. Need 4: droplet effective radius at cloud top and bottom,',...
        'the independent variable, a string describing which independent variable is used,',...
        ' and thermodynamic constraint.', newline])
end

% Check to make sure re values are greater than 0

if any(re_topBottom<0)

    error([newline,'Effective Radius must be greater than 0.', newline])
end


% Check to make sure the thermodynamic constraint is one 1 four
% possibilites

if strcmp(constraint, 'adiabatic')==false && strcmp(constraint, 'subadiabatic_aloft')==false && strcmp(constraint, 'linear_with_z')==false && strcmp(constraint, 'linear_with_tau')==false

    error([newline,'I dont recognize the constraint.', newline])
end


% Check to make sure the independent variable is one of two possibilities

if strcmp(independentVariable, 'altitude')==false && strcmp(independentVariable, 'optical_depth')==false

    error([newline,'I dont recognize the independent variable.', newline])
end




%%



% The boundary values are altered using a powerlaw to get the constraint we
% want




% For the retrieval of r_top, r_bot
if numel(re_topBottom)==2


    % boundary conditions for r as a function of tau
    a0 = @(x) re_topBottom(1)^(2 + 3/x);
    a1 = @(x) re_topBottom(1)^(2 + 3/x) - re_topBottom(2)^(2 + 3/x);

    % boundary conditions for r as a function of z

    b0 = @(x) re_topBottom(2)^(3/x);
    b1 = @(x) re_topBottom(1)^(3/x) - re_topBottom(2)^(3/x);

    if strcmp(constraint,'subadiabatic_aloft')

        % if the profile chosen is subadiabatic aloft, then 0<x<1
        x = 1/2;


        if strcmp(independentVariable,'optical_depth')==true
            re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));

        elseif strcmp(independentVariable,'altitude')==true

            h = max(zT) - min(zT);  % thickness
            h0 = min(zT);           % base height

            re = (b0(x) + b1(x) * (zT - h0)./h).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud

        end


    elseif strcmp(constraint,'adiabatic')

        % if the profile chosen is adiabatic, then x=1
        x = 1;


        if strcmp(independentVariable,'optical_depth')==true

            re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));

        elseif strcmp(independentVariable,'altitude')==true

            h = max(zT) - min(zT);  % thickness
            h0 = min(zT);           % base height

            re = (b0(x) + b1(x) * (zT - h0)./h).^(x/3);              % droplet profile in geometric coordniate system for an adiabatic cloud
        end


    elseif strcmp(constraint,'linear_with_z')

        % linear with z needs an x of 3
        x = 3;

        if strcmp(independentVariable,'optical_depth')==true
            re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));

        elseif strcmp(independentVariable,'altitude')==true

            h = max(zT) - min(zT);
            h0 = min(zT);           % base height

            re = (b0(x) + b1(x) * (zT - h0)./h).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
        end


    elseif strcmp(constraint,'linear_with_tau')

        % linear with tau needs an x of -3
        x = -3;

        if strcmp(independentVariable,'optical_depth')==true
            re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));

        elseif strcmp(independentVariable,'altitude')==true

            h = max(zT) - min(zT);
            h0 = min(zT);           % base height

            re = (b0(x) + b1(x) * (zT - h0)./h).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
        end


    else

        error('I dont recognize the droplet profile you want!')

    end






% For the retrieval of r_top, r_middle, r_bot
elseif numel(re_topBottom)==3



    if strcmp(constraint,'linear_with_tau')

        % We want to allow for two different slopes bewteen r_bot and r_middle
        % and r_middle and r_top


        % boundary conditions for r as a function of tau
        a0_baseMiddle = re_topBottom(2);
        a0_middleTop = re_topBottom(1);

        a1_baseMiddle = re_topBottom(2) - re_topBottom(3);
        a1_middleTop = re_topBottom(1) - re_topBottom(2);


        % boundary conditions for r as a function of z

        b0_baseMiddle = @(x) re_topBottom(3)^(3/x);
        b1_baseMiddle = @(x) re_topBottom(2)^(3/x) - re_topBottom(3)^(3/x);

        b0_middleTop = @(x) re_topBottom(2)^(3/x);
        b1_middleTop = @(x) re_topBottom(1)^(3/x) - re_topBottom(2)^(3/x);



        % linear with tau needs an x of -3
        x = -3;

        if strcmp(independentVariable,'optical_depth')==true


            tau_c = max(zT);
            half_tauC = tau_c/2;
            [~, idx_half_tauC] = min(abs(zT - half_tauC));

            % This computes the droplet profile from cloud base to cloud
            % middle
            re_middleToBase = a0_baseMiddle - a1_baseMiddle * (zT(1:idx_half_tauC)./zT(idx_half_tauC));

            re_topToMiddle = a0_middleTop - a1_middleTop * ((zT(idx_half_tauC:end) - zT(idx_half_tauC))./(zT(end)-zT(idx_half_tauC)));


            % combine these! And don't forgot to skip the repeated value at
            % cloud middle
            % ** This vector starts at cloud top (largest value)
            re = [re_topToMiddle, re_middleToBase(2:end)]; % microns


        elseif strcmp(independentVariable,'altitude')==true

            % let's break this up into two segments
            % first, from base to middle
            h_top = max(zT);            % cloud top height
            h_middle = (max(zT) - min(zT))/2 + min(zT);          % middle height
            h_base = min(zT);           % base height

            [~, idx_middle] = min(abs(zT - h_middle));

            h_baseToMiddle = zT(idx_middle) - h_base;  % cloud thickness from base to middle
            

            h_middleToTop = h_top - zT(idx_middle);  % cloud thickness from middle to top
            

            re_baseToMiddle = (b0_baseMiddle(x) + b1_baseMiddle(x) * (zT(1:idx_middle) - h_base)./h_baseToMiddle).^(x/3);                % droplet profile in geometric coordniate system for an adiabatic cloud

            re_middleToTop = (b0_middleTop(x) + b1_middleTop(x) * (zT(idx_middle:end) - zT(idx_middle))./h_middleToTop).^(x/3);                % droplet profile in geometric coordniate system for an adiabatic cloud


            % combine these! And don't forgot to skip the repeated value at
            % cloud middle
            re = [re_baseToMiddle, re_middleToTop(2:end)]; % microns
        end




    end

end



% define as column vector
re = re';




end
