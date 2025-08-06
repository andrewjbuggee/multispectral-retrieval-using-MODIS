%% ----- READ IN L2 MODIS CLOUD DATA -----

% this function will read in L2 MODIS data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf


% ---- Description of the different fields within the HDF File -----




% By Andrew J. Buggee
%%

function [cloud] = readMODIS_L2_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retreive the scales and offsets for effective particle radius and optical
% thickness

cloudProp_info = hdfinfo(fileName);


% -----------------------------------------------------
% -------- PULL SCALE FACTORS AND OFFSETS -------------
% -----------------------------------------------------

% extract effective radius info first
% These are reported in units of microns!
effectiveRadius17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(5).Value; %  
effectiveRadius17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(6).Value;

% effective radius uncertainty for bands 1 and 7
% These values are reported as percentages of the retireval!
effectiveRadius_uncertainty_17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(5).Value;
effectiveRadius_uncertainty_17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(6).Value;

% Effective radius scale factor and offset for bands 1 and 6
effectiveRadius16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(5).Value; %  
effectiveRadius16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(6).Value;


% effective radius uncertainty for bands 1 and 6
% uncertainties are listed as percents
effectiveRadius_uncertainty_16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(5).Value;
effectiveRadius_uncertainty_16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(6).Value;



% extract the optical thickness scales and offset
optThickness17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(5).Value; % 
optThickness17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(6).Value;


% optical thickness uncertainty for bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(94).Attributes(5).Value;
optThickness_uncertainty_17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(94).Attributes(6).Value;


% optical thickness scale factors and offsets for bands 1 and 6
optThickness16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(5).Value; %  output will be a cell array
optThickness16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(6).Value;

% optical thickness uncertainty for bands 1 and 6
% uncertainties are listed as percents
optThickness_uncertainty_16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(95).Attributes(5).Value;
optThickness_uncertainty_16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(95).Attributes(6).Value;

% extract the cloud top height at 1km resolution scales and offset
cloudTopHeight_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(58).Attributes(5).Value; %  output will be a cell array
cloudTopHeight_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(58).Attributes(6).Value;

% extract the cloud top pressure at 1km resolution scales and offset
cloudTopPressure_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(57).Attributes(5).Value; %  output will be a cell array
cloudTopPressure_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(57).Attributes(6).Value;

% extract the cloud top temperature at 1km resolution scales and offset
cloudTopTemperature_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(59).Attributes(5).Value; %  output will be a cell array
cloudTopTemperature_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(59).Attributes(6).Value;

% extract the liquid water path at 1km resolution scales and offset
columnWaterPath_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(83).Attributes(5).Value; %  output will be a cell array
columnWaterPath_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(83).Attributes(6).Value;

% extract the cloud phase used in Optical Thickness/Effective Radius determination -  scales and offset
% The values in this SDS are set to mean the following:                              
% 0 -- cloud mask undetermined                                                       
% 1 -- clear sky                                                                     
% 2 -- liquid water cloud                                                            
% 3 -- ice cloud                                                                     
% 4 -- undetermined phase cloud (but retrieval is attempted as  liquid water)        

cloudPhase_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(105).Attributes(5).Value; %  output will be a cell array
cloudPhase_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(105).Attributes(6).Value;


% extract the cloud fraction
cloudFraction_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(35).Attributes(5).Value; %  output will be a cell array
cloudFraction_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(35).Attributes(6).Value;


% extract the cloud fraction - day time only
cloudFractionDay_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(39).Attributes(5).Value; %  output will be a cell array
cloudFractionDay_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(39).Attributes(6).Value;

% extract the above cloud water vapor
% determined by the 0.94 micron channel
% only applies to pixels above ocean with an optical thickness greater>5
% units are in cm
aboveCloudWaterVapor_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(103).Attributes(5).Value; %  output will be a cell array
aboveCloudWaterVapor_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(103).Attributes(6).Value;




% -----------------------------------------------------
% ---------------- PULL THE DATA SET ------------------
% -----------------------------------------------------

% uncertainties are listed as percents


% extract the effective radius data using bands 1 and 7
effectiveRadius_17 = hdfread(fileName,'Cloud_Effective_Radius');

% extract the effective radius uncertainty bands 1 and 7
% uncertainties are listed as percents
effectiveRadius_uncertainty_17 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty');

% extract the effective radius data using bands 1 and 6
effectiveRadius_16 = hdfread(fileName,'Cloud_Effective_Radius_16');

% extract the effective radius uncertainty bands 1 and 6
% uncertainties are listed as percents
effectiveRadius_uncertainty_16 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty_16');


% extract the Optical thickness data using bands 1 and 7
opticalThickness_17 = hdfread(fileName,'Cloud_Optical_Thickness');

% extract the optical thickness uncertainty bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_17 = hdfread(fileName,'Cloud_Optical_Thickness_Uncertainty');



% extract the Optical thickness data using bands 1 and 6
opticalThickness_16 = hdfread(fileName,'Cloud_Optical_Thickness_16');

% extract the optical thickness uncertainty bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_16 = hdfread(fileName,'Cloud_Optical_Thickness_Uncertainty_16');


% extract the cloud top geopotential height at 1km resolution - units (meters)
cloudTop_geopotentialHeight = hdfread(fileName,'cloud_top_height_1km');

% extract the cloud phase - 
% The values are:
%   0 - no phase result
%   1 - no phase result
%   2 - liquid water
%   3 - ice
%   4 - undetermined

% extract the cloud top pressure - units (hpa)
cloudTop_pressure = hdfread(fileName,'cloud_top_pressure_1km');

% extract the cloud top temperature - units (k)
cloudTop_temperature = hdfread(fileName,'cloud_top_temperature_1km');

% extract the cloud water height for band 7 and either 1 or 2 or 5
columnWaterPath = hdfread(fileName, 'Cloud_Water_Path');            % g/m^2

% extract the cloud phase - 
% The values are:
%   0 - no phase result
%   1 - no phase result
%   2 - liquid water
%   3 - ice
%   4 - undetermined

cloudPhase = hdfread(fileName,'Cloud_Phase_Optical_Properties');


% Read the cloud_mask_SPI variable which contains the sub-pixel
% heterogeneity index. Larger values have been shown the have more
% retrieval bias and increased rates of retreival failure

subPix_heteroIndex = hdfread(fileName, 'Cloud_Mask_SPI');

subPix_heteroIndex_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(118).Attributes(5).Value;
subPix_heterIndex_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(118).Attributes(6).Value;


% extract the cloud fraction data
% values are listed as a percent between 0 and 100
cloudFraction = hdfread(fileName, 'Cloud_Fraction');

% extract the cloud fraction data for day time only
% values are listed as a percent between 0 and 100
cloudFractionDay = hdfread(fileName, 'Cloud_Fraction_Day');


% extract the above cloud water vapor
% values in units of cm - the column water vapor if condensed to a square
% cm, this is the height of the column of condensed water
aboveCloudWaterVapor = hdfread(fileName, 'Above_Cloud_Water_Vapor_094');





%% --- Convert Data ---


cloud.effRadius17 = scalesOffsets2Matrix(effectiveRadius_17,effectiveRadius17_scales,effectiveRadius17_offset);
cloud.effRad_uncert_17 = scalesOffsets2Matrix(effectiveRadius_uncertainty_17,effectiveRadius_uncertainty_17_scales,effectiveRadius_uncertainty_17_offset);
cloud.optThickness17 = scalesOffsets2Matrix(opticalThickness_17,optThickness17_scales,optThickness17_offset);
cloud.optThickness_uncert_17 = scalesOffsets2Matrix(optThickness_uncertainty_17,optThickness_uncertainty_17_scales,optThickness_uncertainty_17_offset);

cloud.effRadius16 = scalesOffsets2Matrix(effectiveRadius_16,effectiveRadius16_scales,effectiveRadius16_offset);
cloud.effRad_uncert_16 = scalesOffsets2Matrix(effectiveRadius_uncertainty_16,effectiveRadius_uncertainty_16_scales,effectiveRadius_uncertainty_16_offset);
cloud.optThickness16 = scalesOffsets2Matrix(opticalThickness_16,optThickness16_scales,optThickness16_offset);
cloud.optThickness_uncert_16 = scalesOffsets2Matrix(optThickness_uncertainty_16,optThickness_uncertainty_16_scales,optThickness_uncertainty_16_offset);

cloud.topHeight = scalesOffsets2Matrix(cloudTop_geopotentialHeight,cloudTopHeight_scales,cloudTopHeight_offset);
cloud.topPressure = scalesOffsets2Matrix(cloudTop_pressure,cloudTopPressure_scales,cloudTopPressure_offset);
cloud.topTemperature = scalesOffsets2Matrix(cloudTop_temperature,cloudTopTemperature_scales,cloudTopTemperature_offset);

cloud.phase = scalesOffsets2Matrix(cloudPhase,cloudPhase_scales,cloudPhase_offset);

cloud.lwp = scalesOffsets2Matrix(columnWaterPath, columnWaterPath_scales, columnWaterPath_offset);

cloud.SPI = scalesOffsets2Matrix(subPix_heteroIndex,subPix_heteroIndex_scales,subPix_heterIndex_offset);

cloud.fraction = scalesOffsets2Matrix(cloudFraction,cloudFraction_scales,cloudFraction_offset);

cloud.fraction_day = scalesOffsets2Matrix(cloudFractionDay,cloudFractionDay_scales,cloudFractionDay_offset);


cloud.aboveWaterVaporCol = scalesOffsets2Matrix(aboveCloudWaterVapor,aboveCloudWaterVapor_scales,aboveCloudWaterVapor_offset);


%% ---- Convert null results to nan's -----

% Cloud effective radius cannot be less than 0. Occasionally MODIS data
% will contain negative values for measurements that need to be thrown out.
% Convert these to nans

cloud.effRadius17(cloud.effRadius17<-99) = nan;
cloud.effRadius16(cloud.effRadius16<-99) = nan;

% Cloud optical thickness cannot be less than 0. Occasionally MODIS data
% will contain negative values for measurements that need to be thrown out.
% Convert these to nans

cloud.optThickness17(cloud.optThickness17<-99) = nan;
cloud.optThickness16(cloud.optThickness16<-99) = nan; 

% Above cloud water vapor is measured in cm and can only be greater than 0.
% All values less than 0 represent a retrieval for a cloud with an optical
% thickness less than 5 or a pixel over ocean

cloud.aboveWaterVaporCol(cloud.aboveWaterVaporCol<0) = nan;





end





