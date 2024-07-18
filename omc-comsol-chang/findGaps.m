% This function finds any complete band gaps in a given band structure and
% returns a list of their mid-gap frequencies and size of the bandgaps. The
% input structure ds is assumed to have four fields: EzEy, EzOy, OzEy, and
% OzOy, corresponding to the different mode symmetries. Each field is
% assumed to have a matrix F which stores each band in a different column.
%
% Sean Meenehan, 10/7/11

function [midGap,gapSize] = findGaps(ds)

% Assemble a matrix of all bands, then find the max and min of each band
bands = ds.F;
sortBands = sort(bands);
bMin = sortBands(1,:);
bMax = sortBands(end,:);

midGap = [];
gapSize = [];
for k = 1:length(bMin)
    % To test for a band gap, we look at the minimum of each band and
    % see whether there exist any bands for which it lies between the minimum
    % and maximum. If all other bands are either completely above or
    % completely below the minimum, then it must be the top of a complete
    % band gap
    currMin = bMin(k);
    bandTest = 1; 
    for t = 1:length(bMin)
        if t ~= k
            testMin = bMin(t);
            testMax = bMax(t);
            if (currMin > testMin && currMin < testMax)
                bandTest = 0;
            end
        end
    end
    
    if bandTest
        maxInds = find(bMax < currMin);
        if ~isempty(maxInds)
            maxVal = max(bMax(maxInds));
            midGap(end+1) = (currMin+maxVal)/2;
            gapSize(end+1) = currMin-maxVal;
        end
    end
end
