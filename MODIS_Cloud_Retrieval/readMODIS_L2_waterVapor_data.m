%% ----- READ IN L2 MODIS WATER VAPOR DATA -----

% this function will read in L2 MODIS water vapor data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf

% This data is designated with the tag '-05' such as 'MOD05_L2...' for the
% terra platform and 'MYD05_L2...' for the Aqua platform.


% ---- Description of the different fields within the HDF File -----




% By Andrew J. Buggee
%%

function [vapor] = readMODIS_L2_waterVapor_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retreive the scales and offsets for effective particle radius and optical
% thickness

waterVapor_info = hdfinfo(fileName);


% -----------------------------------------------------
% -------- PULL SCALE FACTORS AND OFFSETS -------------
% -----------------------------------------------------

% extract the near infrared total column water vapor scale and offset
columnWaterVapor_NIR_scales = waterVapor_info.Vgroup.Vgroup(2).SDS(7).Attributes(3).Value; %  
columnWaterVapor_NIR_offset = waterVapor_info.Vgroup.Vgroup(2).SDS(7).Attributes(4).Value;



% extract the near infrared total column water vapor correction factor scale and offset
columnWaterVaporCorrection_NIR_scales = waterVapor_info.Vgroup.Vgroup(2).SDS(8).Attributes(3).Value; %  
columnWaterVaporCorrection_NIR_offset = waterVapor_info.Vgroup.Vgroup(2).SDS(8).Attributes(4).Value;



% extract the infrared total column water vapor scale and offset
columnWaterVapor_IR_scales = waterVapor_info.Vgroup.Vgroup(2).SDS(9).Attributes(5).Value; %  
columnWaterVapor_IR_offset = waterVapor_info.Vgroup.Vgroup(2).SDS(9).Attributes(6).Value;






% -----------------------------------------------------
% ---------------- PULL THE DATA SET ------------------
% -----------------------------------------------------

% uncertainties are listed as percents


% extract the near IR water vapor retrieval
column_waterVapor_NIR = hdfread(fileName,'Water_Vapor_Near_Infrared');      % cm


% extract the near IR water vapor retrieval correction
column_waterVapor_NIR_correction = hdfread(fileName,'Water_Vapor_Correction_Factors');  

% extract the IR water vapor retrieval
column_waterVapor_IR = hdfread(fileName,'Water_Vapor_Infrared');      % cm



%% --- Convert Data ---


vapor.col_nir = scalesOffsets2Matrix(column_waterVapor_NIR,columnWaterVapor_NIR_scales,columnWaterVapor_NIR_offset);        % cm

vapor.col_nir_correction = scalesOffsets2Matrix(column_waterVapor_NIR_correction,columnWaterVaporCorrection_NIR_scales,columnWaterVaporCorrection_NIR_offset); 
    
vapor.col_ir = scalesOffsets2Matrix(column_waterVapor_IR,columnWaterVapor_IR_scales,columnWaterVapor_IR_offset);         % cm





%% ---- Convert null results to nan's -----

% Cloud effective radius cannot be less than 0. Occasionally MODIS data
% will contain negative values for measurements that need to be thrown out.
% Convert these to nans

vapor.col_nir(vapor.col_nir<-99) = nan;
vapor.col_ir(vapor.col_ir<-99) = nan;







end





