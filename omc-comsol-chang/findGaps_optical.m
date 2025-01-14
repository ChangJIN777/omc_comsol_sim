% This function finds any complete band gaps in a given band structure and
% returns a list of their mid-gap frequencies and size of the bandgaps. The
% input structure ds is assumed to have four fields: EzEy, EzOy, OzEy, and
% OzOy, corresponding to the different mode symmetries. Each field is
% assumed to have a matrix F which stores each band in a different column.
%
% Sean Meenehan, 10/7/11

function [midGap,gapSize] = findGaps_optical(ds)

% Assemble a matrix of all bands, then find the max and min of each band
% each col corresponds to each band
bands = ds.F;
sortedBands = sort(bands); % each col is now sorted from min to max freq, then NaNs
bMin = min(sortedBands); %find min and max of each col
bMax = max(sortedBands);
bMinTrueIdx = find(~isnan(bMin)); % find bands where both min and max of band are not NaN
bMaxTrueIdx = find(~isnan(bMax));
bTrueIdx = intersect(bMinTrueIdx,bMaxTrueIdx);
bMin = bMin(bTrueIdx);
bMax = bMax(bTrueIdx);

% 
% bMinTmp = sortBands(1,:);
% bMaxTmp = sortBands(end,:);
% 
% %find min and max from each band excluding NaN
% actualbMin = bMinTmp(~isnan(bMinTmp));
% if length(actualbMin) ~= length(bMinTmp)
%     bMin = zeros(1,length(actualbMin));
%     bMax = zeros(1,length(actualbMin));
%     for bi = 1:length(actualbMin)
%         band = bands(:,bi);
%         actualBand = band(~isnan(band));
%         bMin(bi) = min(actualBand);
%         bMax(bi) = max(actualBand);
%     end
% else
%     bMin = bMinTmp;
%     bMax = bMaxTmp;
% end

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