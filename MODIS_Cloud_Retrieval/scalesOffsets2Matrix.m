%% ----- Convert the Scales and Offsets into Matrices -----

% Convert the scales and offset vectors into the appropriate sized matrix
% for easy conversion between scaled integer values and other quantities of
% interest. Then we comput the new matrix, which is either the radiance, or
% the reflectance. Both conversions have the same form:

%       new(T,FS,B) = scales(B) * (scaled_integers(T,FS,B) - offsets(B)

% where T = track, FS = frame and sample, and B = band

% By Andrew J. Buggee
%%

function [radOrRefMatrix] = scalesOffsets2Matrix(scaledIntegerMatrix,scales,offsets)

% convert everything to type double

numRows = size(scaledIntegerMatrix,1);
numCols = size(scaledIntegerMatrix,2);
numBands = size(scaledIntegerMatrix,3);

if length(offsets)==1
    
    if length(scales)>1
        
        scalesMat = repmat(reshape(scales,1,1,[]),numRows,numCols);
        offsetsMat = repmat(offsets,numRows,numCols);
        
    elseif length(scales)==1
        
        scalesMat = repmat(scales,numRows,numCols);
        offsetsMat = repmat(offsets,numRows,numCols);
        
    end
    
else
    
    if length(scales)>1
        
        scalesMat = repmat(reshape(scales,1,1,[]),numRows,numCols);
        offsetsMat = repmat(reshape(offsets,1,1,[]),numRows,numCols);
        
    elseif length(scales)==1
        
        scalesMat = repmat(scales,numRows,numCols);
        offsetsMat = repmat(reshape(offsets,1,1,[]),numRows,numCols);
        
    end
    
    
    
    
    
end

radOrRefMatrix = double(scalesMat) .* (double(scaledIntegerMatrix) - double(offsetsMat));

end



