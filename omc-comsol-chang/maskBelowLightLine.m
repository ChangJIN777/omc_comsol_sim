function [Fmasked,below] = maskBelowLightLine(F,lightline,tol)
%MASKBELOWLIGHTLINE Keep only modes strictly below the light line, with margin.
%
% Returns the band matrix with every radiative point replaced by NaN, plus the
% logical mask itself. Both solveOpticalBands and test_plotOpticalBands call
% this so the set of points feeding the gap search cannot differ between a live
% run and a re-analysis of its saved data.
%
% INPUTS
%   F          [nk x nbands] frequencies in Hz.
%   lightline  [nk x 1] light line in Hz at each k-point, from opticalLightLine.
%   tol        Relative margin, dimensionless. A point counts as guided only if
%              f < lightline*(1-tol). Optional, default 1e-3 (0.1%).
%
% OUTPUTS
%   Fmasked    F with non-guided entries set to NaN.
%   below      [nk x nbands] logical, true where the point is kept.
%
% WHY A MARGIN AND NOT JUST '<'
% The bare test f < lightline already excludes a point exactly ON the line, but
% exact equality between two computed doubles essentially never occurs, so that
% test does nothing to points sitting on the line to within solver noise. Those
% are the least trustworthy points in the whole diagram: a mode at the light
% line is only marginally bound, its classification flips with mesh refinement,
% and because findGaps_belowLightLine works on the SET of surviving frequencies,
% one such point sitting inside an otherwise-empty region destroys the gap it
% falls in. The margin deactivates that band of ambiguity explicitly rather than
% leaving it to rounding.
%
% Set tol = 0 to recover the previous strict-inequality behaviour exactly.
%
% NOTE ON THE OLD IDIOM: the callers previously did
%     TE.F = F.*below;  TE.F(TE.F==0) = NaN;
% which uses 0 as a sentinel for "masked" while 0 is also a legitimate
% frequency - a genuine 0 Hz mode, at Gamma say, was silently discarded as if it
% were radiative. Assigning NaN directly, as here, keeps the mask and the data
% separate.
%
% EXAMPLE
%   ll = opticalLightLine(OpticalBand,P);
%   [Fm,below] = maskBelowLightLine(OpticalBand.F,ll,1e-3);
%   [midGap,gapSize] = findGaps_belowLightLine(Fm,20e12,100e12);

narginchk(2,3);
if nargin < 3 || isempty(tol)
    tol = 1e-3;
end
validateattributes(F,{'numeric'},{'2d'},mfilename,'F');
validateattributes(tol,{'numeric'},{'scalar','nonnegative','finite','<',1}, ...
    mfilename,'tol');

% lightline is [nk x 1] and F is [nk x nbands]; implicit expansion compares each
% band against its own k-point's light line. Column orientation is forced so a
% row-vector lightline cannot silently broadcast along the wrong dimension.
lightline = lightline(:);
if numel(lightline) ~= size(F,1)
    error('maskBelowLightLine:sizeMismatch', ...
        ['lightline has %d entries but F has %d k-point rows. They must ' ...
         'match one-to-one.'],numel(lightline),size(F,1));
end

below          = F < lightline*(1-tol);
Fmasked        = F;
Fmasked(~below) = NaN;
end
