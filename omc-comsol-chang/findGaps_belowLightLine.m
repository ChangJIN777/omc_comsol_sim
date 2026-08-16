function [midGap,gapSize,gapEdges,maxRejected] = findGaps_belowLightLine(F,minGap)
%FINDGAPS_BELOWLIGHTLINE Frequency regions containing no guided data point.
%
% Treats the band structure as an unordered CLOUD of frequencies and reports
% every interval between consecutive occupied frequencies that is at least
% minGap wide. A reported gap is, by construction, a frequency band that no
% surviving data point falls inside.
%
% This replaces findGaps_optical on the optical path. findGaps_optical is left
% untouched on disk but now has NO callers - the mechanical pipeline uses
% findGaps, a different function. The two answer different questions:
%
%   findGaps_optical         collapses each BAND to its [min,max] interval and
%                            looks for band minima that no other band's interval
%                            straddles. Correct only if every column of F is a
%                            consistently-ordered level; scrambled columns widen
%                            the intervals until every band straddles every
%                            other and all gaps silently disappear.
%
%   findGaps_belowLightLine  ignores columns entirely. Only the SET of occupied
%                            frequencies matters, so the result cannot be
%                            affected by how points are distributed among bands.
%                            That matters here: runOpticalBand_* assigns a point
%                            to a band purely by its position in COMSOL's
%                            eigenvalue list at that k-point, with no mode
%                            tracking and no sort, so column identity is not
%                            something to lean on.
%
% INPUTS
%   F        [nk x nbands] frequencies in Hz. NaN marks a point that must not
%            count as occupied - in this project, a mode above the light line,
%            masked out by solveOpticalBands before the call. Inf and complex
%            entries are treated the same way as NaN.
%   minGap   Minimum width in Hz for a region to be reported. Optional,
%            default 0 (report every gap, however narrow).
%
% OUTPUTS
%   midGap       [n x 1] centre frequency of each region, Hz.
%   gapSize      [n x 1] width of each region, Hz.
%   gapEdges     [n x 2] lower and upper edge of each region, Hz. The upper
%                edge is the lowest guided point above the gap, the lower edge
%                the highest guided point below it.
%   maxRejected  Width of the widest region that failed minGap, Hz, or [] if
%                none did. Reported so "no gaps" can be distinguished from
%                "gaps, but all narrower than the threshold" without re-running.
%
% All outputs are column vectors sorted by ascending frequency, and are empty
% (0x1, or 0x2 for gapEdges) when nothing qualifies.
%
% WHAT IS DELIBERATELY NOT REPORTED
%   - The region below the lowest guided point. That is not a bandgap, it is
%     the space under the fundamental band.
%   - The region above the highest guided point. That is an artefact of solving
%     a finite number of bands (P.nbands): more bands exist up there, they were
%     simply never computed, so the apparent void is unbounded and meaningless.
%
% LIMITATION worth knowing before quoting a number: gaps are measured from the
% SAMPLED k-points only. A band extremum falling between two sampled k-points is
% invisible, so a reported gap is an upper bound on the true one and can be
% spurious entirely if the sampling is coarse. With P.kpts = 9 the circuit is
% sampled at 27 points. Increase kpts before trusting a marginal gap.
%
% EXAMPLE
%   TEbelow = OpticalBand.F < lightline;
%   Fmasked = OpticalBand.F;  Fmasked(~TEbelow) = NaN;
%   [midGap,gapSize] = findGaps_belowLightLine(Fmasked,10e12);
%   fprintf('%.2f THz gap at %.2f THz\n',[gapSize gapMid].'*1e-12);

narginchk(1,2);
if nargin < 2 || isempty(minGap)
    minGap = 0;
end
validateattributes(F,{'numeric'},{'2d'},mfilename,'F');
validateattributes(minGap,{'numeric'},{'scalar','nonnegative','finite'}, ...
    mfilename,'minGap');

midGap      = zeros(0,1);
gapSize     = zeros(0,1);
gapEdges    = zeros(0,2);
maxRejected = [];

% Flatten to the set of occupied frequencies. isfinite drops NaN (above the
% light line) and Inf; real() guards against a complex eigenfrequency leaking
% through, which would otherwise make sort and diff behave by magnitude.
f = real(F(:));
f = f(isfinite(f));
if numel(f) < 2
    return          % nothing, or a single point: no interval to speak of
end

f = sort(f);
d = diff(f);        % width of every interval between consecutive occupied freqs

qualifies = d >= minGap;
if ~any(qualifies)
    if ~isempty(d)
        maxRejected = max(d);
    end
    return
end

idx      = find(qualifies);
gapEdges = [f(idx), f(idx+1)];
gapSize  = d(idx);
midGap   = (f(idx) + f(idx+1))/2;

rejected = d(~qualifies);
if ~isempty(rejected)
    maxRejected = max(rejected);
end
end
