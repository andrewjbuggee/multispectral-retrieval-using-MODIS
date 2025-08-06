% Convert the MODIS data from format [band,acrossTrack,alongTrack]] to
% [acrossTrack,AlongTrack,band]. This makes the data easier to work with in
% Matlab

function [newData] = bandAcrossAlong2AcrossAlongBand(data)

% split the 250 meter earth view data in (row,column,band) format

numRows = size(data,2); % across track direction, taken using scan mirror. calculated by num_scans * num_pixels
numCols = size(data,3); % along track direction, taken using num_frames*num_samples
numBands = size(data,1); % bands

newData = zeros(numRows,numCols,numBands);

for ii = 1:numBands
    
    newData(:,:,ii) = reshape(data(ii,:,:),numRows,numCols);
    
end

end

