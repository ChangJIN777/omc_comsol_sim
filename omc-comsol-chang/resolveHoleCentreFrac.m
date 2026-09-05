function f = resolveHoleCentreFrac(P)
%RESOLVEHOLECENTREFRAC Where the boomerang hole sits inside the rhombic cell.
%
% Returns the fraction f such that the hole centre is
%
%     c = f*(a1 + a2)  =  f*[3*a/2, a*sqrt(3)/2]
%
% with a1 = (a,0) and a2 = (a/2, a*sqrt(3)/2) the primitive vectors of the
% hexagonal lattice. f = 0.5 puts the hole on the cell centroid; smaller f
% slides it toward the origin along the cell diagonal, which is the 210 deg arm
% direction.
%
% DEFAULT f = 0.4, not 0.5. The centroid is NOT the position that allows the
% largest arm length r before the hole crosses the cell boundary, because the
% hole is 3-fold symmetric (arms at 90/210/330 deg) while the rhombic cell is
% only 2-fold - the constraints are unbalanced. Sliding down the diagonal trades
% clearance from the two down-pointing arms, which have room to spare, to the
% single up-pointing arm, which is the one that binds. Measured for a = 730 nm:
%
%     f      max r before the hole leaves the cell
%     0.50    316.1 nm   (centroid)
%     0.40    ~379 nm
%     0.39    384.9 nm   (optimum at w = 125 nm)
%
% The optimum drifts with w - 0.367 at w = 75 nm, 0.425 at w = 200 nm - so 0.4
% is the round compromise across a realistic sweep rather than the exact optimum
% for any single w. Holes in neighbouring cells merge at r ~ 400 nm regardless
% of f, so f = 0.4 recovers nearly all the range physically available.
%
% WHAT THIS DOES AND DOES NOT CHANGE
%   Does NOT change the physics. Translating the motif inside a periodic cell
%   produces an identical crystal - same lattice, same motif, only the
%   bookkeeping moves. Band structures are unaffected.
%
%   DOES change how large r can be before the drawn hole crosses the cell
%   boundary. Past that point Compose clips the hole and the 'inside'
%   BoxSelections in buildBoomerangUnitCell start swallowing interior faces,
%   putting Floquet conditions on the wrong surfaces - see the warning at
%   buildBoomerangUnitCell.m:251-255. That is the failure this pushes back.
%
% Set P.holeCentreFrac = 0.5 to restore the original centred geometry.
%
% INPUT
%   P   parameter struct; reads the optional field holeCentreFrac.
%
% OUTPUT
%   f   scalar in (0,1).
%
% Shared by buildBoomerangUnitCell (what COMSOL builds) and calcFillingFactor
% (what the pre-solve gate measures) so the two cannot disagree about where the
% hole is. Adding a third consumer means calling this, not re-declaring 0.4.

narginchk(1,1);

if isfield(P,'holeCentreFrac') && ~isempty(P.holeCentreFrac)
    f = P.holeCentreFrac;
else
    f = 0.4;
end

% Bounds are open: f = 0 or 1 puts the hole centre on a cell vertex, where every
% arm immediately leaves the cell. Caught here rather than as an opaque COMSOL
% geometry failure several calls deeper.
validateattributes(f,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<',1}, mfilename, 'P.holeCentreFrac');
end
