%% ----- Write output file names for uvspec -----




% By Andrew J. Buggee

%%


function [outputNames] = writeOutputNames(inputNames)

numPixels = size(inputNames,1);
num_radii = size(inputNames,2);
num_tau = size(inputNames,3);
num_bands = size(inputNames,4);

outputNames = cell(size(inputNames));

for pp = 1:numPixels
    
    for rr = 1:num_radii
        
        for tt = 1:num_tau
            
            for bb = 1:num_bands
                
                outputNames{pp,rr,tt,bb} = ['OUTPUT_',inputNames{pp,rr,tt,bb}(1:end-4)];
                
                
            end
            
        end
        
    end
    
end



end