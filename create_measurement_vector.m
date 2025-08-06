% --- Create the measurement vector needed for Rodgers solution of gaussian
% model pdf and measurement pdf ---


% By Andrew J. Buggee

%%

function y = create_measurement_vector(modis, GN_inputs, pixels2use)


% We will retireve cloud properties for n pixels, designated by the
% Gauss_Newton input structure
num_pixels = GN_inputs.numPixels2Calculate;

% Use the bands specified by GN_inputs.bands2use.  Most of the time this is
% the first 7 bands

y = zeros(length(GN_inputs.bands2use),num_pixels);

for pp = 1:num_pixels
    
    
    row = pixels2use.res1km.row(pp);
    col = pixels2use.res1km.col(pp);
    
    y(:,pp) = reshape(modis.EV1km.reflectance(row,col,GN_inputs.bands2use),[],1);

    % ----------------------------------------------------------------------
    % The data from 11-11-2008 at 18:50 measured erroneous values in the 1.6
    % micron channel. If using this data, lets ignore this measurement. This
    % may happen with other measurements. Always check the reflectance
    % uncertainty. If this value is above 10%, let's issue a warning and remove
    % it from use in the retrieval. Afterall, uncertainty is the MOST important
    % parameter when retrieval droplet profiles
    
    idx_above10 = reshape(modis.EV1km.reflectanceUncert(row, col, GN_inputs.bands2use), [],1) > 10;     % uncertainty listed as a percent

    % This will throw an error if you try to use more than 1 pixel at a
    % time. For now, rather than change my algorithm to handle cases
    % where some measurements shouldn't be used, just run the offending
    % pixel on it's own, rather than in a loop with other pixels.
    if any(idx_above10)
        y(idx_above10, pp) = nan;
    
        warning([newline, 'Reflectance measured at pixel [r,c] = [', num2str(row), ', ', num2str(col), '], has ',...
            'an uncertainty above 10% at band number ', num2str(find(idx_above10)),'.', newline,...
            'The reflectance at this band has been set to NaN.', newline])
    end
    % ----------------------------------------------------------------------
    
end







end