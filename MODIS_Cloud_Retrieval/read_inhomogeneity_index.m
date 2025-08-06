%% ----- READ IN MOD/MYD35 file and extract the sub pixel inhomoegeneity index -----


% By Andrew John Buggee

%%

function [inhomogeneity_index] = read_inhomogeneity_index(fileName)

% read cloud mask and spectral test file
cloudMask_info = hdfinfo(fileName);

% extract the sub-pixel-inhomogeneity index
% First dimension is the inhomogeneity index for band 1 (650 nm)
% Second dimension is the inhomogeneity index for band 2 (860 nm)
% H = std(R(250m))/mean(R(250m))
subPixel_inhomogeneity_index = hdfread(fileName,'Cloud_Mask_SPI');      % percent

% extract the inhomogeneity index scale and offset
inhomogeneity_index_scales = cloudMask_info.Vgroup.Vgroup(2).SDS(6).Attributes(5).Value; %  
inhomogeneity_index_offset = cloudMask_info.Vgroup.Vgroup(2).SDS(6).Attributes(6).Value;


%% Convert data using the scales and offsets

 inhomogeneity_index = scalesOffsets2Matrix(subPixel_inhomogeneity_index, inhomogeneity_index_scales, inhomogeneity_index_offset);         % cm




end
