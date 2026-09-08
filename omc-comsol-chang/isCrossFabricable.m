function [ok, margin, why] = isCrossFabricable(a, h, w, th, minFeature, r1, r2)
%ISCROSSFABRICABLE  Lithography feasibility of a cross unit cell, no solve.
%
%   [ok, margin, why] = isCrossFabricable(a, h, w, th, minFeature, r1, r2)
%
%   VECTORIZED over a/h/w/th (scalars broadcast), because bayesopt hands an
%   XConstraintFcn a TABLE of many candidate rows at once and expects a logical
%   column back. The scalar path is just the n = 1 case.
%
%   GEOMETRY (verified against buildCrossUnitCell.m, not from the comments in
%   test_CrossUnitCell.m, which describe a different cell):
%     a   square cell side, the lattice constant in BOTH x and y  (:38, [a a])
%     h   LENGTH of each cross arm                          (:43 [h w], :48 [w h])
%     w   WIDTH of each cross arm
%     th  slab thickness along z                                  (:33, :65)
%     r1  fillet radius 1                                              (:137)
%     r2  fillet radius 2                                              (:147)
%
%   The compose formula is 'r_ucell-r1-r2' (:53), so the cross is SUBTRACTED:
%   the solid is a square slab with a cross-shaped VOID. In the xy plane the
%   void is the union of two crossed bars,
%       horizontal arm :  |x| <= h/2 , |y| <= w/2
%       vertical   arm :  |x| <= w/2 , |y| <= h/2
%   so the limits below are written on the SOLID that survives and on the
%   etched gap width, not on the arms as if they were solid.
%
%   OUTPUTS
%     ok     - logical array, true where the design is fabricable
%     margin - worst binding margin in METRES. Negative means violated, and its
%              magnitude says by how much, so a caller can build a GRADED
%              penalty rather than a flat one.
%     why    - cell array naming the tightest rule ('' where ok)
%
%   WHY CLOSED FORM IS ENOUGH HERE. calcFillingFactor.m has to MEASURE the
%   equivalent quantity for the boomerang cell: it rebuilds the void as a
%   polyshape, densifies its boundary, translates it to the six nearest
%   lattice sites and takes the closest approach, because a boomerang void is
%   a complex, possibly off-centre shape on a RHOMBIC lattice, where the
%   thinnest wall only exists once the lattice is applied. The cross void is
%   two axis-aligned rectangles, centred, on a SQUARE lattice, so every
%   clearance has an exact expression and there is nothing a measurement
%   would find that arithmetic does not. The wall formula below was checked
%   against the tiled polygons anyway; see its comment.
%
%   THESE RULES ARE INFERRED FROM THE GEOMETRY CONSTRUCTION, NOT FROM A PROCESS
%   SPEC. Adjust minFeature, or the rules themselves, to match your actual
%   lithography and etch limits before trusting them to prune real designs.
%
%   Example:
%     ok = isCrossFabricable(160e-9, 140e-9, 50e-9, 120e-9, 10e-9, 10e-9, 10e-9)
%
%   See also BUILDCROSSUNITCELL, BAYESOPT_CROSS, CROSS_OPTIMIZE_SWEEP_DIAMOND.

narginchk(5, 7);
if nargin < 6; r1 = 0; end
if nargin < 7; r2 = 0; end

% Broadcast to a common size so scalar/array mixes work.
%
% Done via implicit expansion rather than by comparing size() vectors. The
% obvious max([size(a); size(h); ...]) fails the moment the inputs differ in
% NUMBER of dimensions -- a 21x21x21 grid gives a 1x3 size vector while a
% scalar gives 1x2, and the vertcat errors. That is not hypothetical: it is
% exactly the shape of the call the fixed-thickness path makes, three gridded
% variables against a scalar th.
probe = a + h + w + th;      % also errors clearly if two arrays are incompatible
sz    = size(probe);
a  = a  .* ones(sz);
h  = h  .* ones(sz);
w  = w  .* ones(sz);
th = th .* ones(sz);
r1 = r1 .* ones(sz);           % each radius is bounded separately below,
r2 = r2 .* ones(sz);           % so max(r1,r2) is no longer sufficient

% Absolute tolerance on every comparison. Without it a design sitting exactly
% ON a limit is rejected by floating-point noise: (160e-9 - 140e-9)/2 evaluates
% to 9.9999999999999969e-9 rather than 1e-8, and so "fails" a 10 nm rule by
% 3e-24 m. 1 fm is twenty orders of magnitude below any real process tolerance,
% so it can only ever absorb representation error.
tol = 1e-15;

% THE NARROWEST SOLID FEATURE IS THE WALL BETWEEN THE VOIDS OF TWO ADJACENT
% CELLS, not the distance from an arm end to the cell boundary. The cell
% boundary is a periodic-BC plane, not a piece of geometry -- the material on
% the far side belongs to the next cell and is continuous with this one. Tile
% the square lattice and an arm end at x = h/2 faces the opposite arm end of
% the neighbour at x = a - h/2, so
%
%     wall = (a - h/2) - h/2 = a - h
%
% Verified numerically against the tiled polygons (closest approach between a
% cell's void and its eight nearest translates): a = 160, h = 140 nm gives a
% measured wall of 20.00 nm, exactly a - h.
%
% This rule previously constrained (a-h)/2, the HALF-distance to the boundary,
% making it 2x stricter than the real process limit and quietly halving the
% usable search box. It never admitted anything unbuildable; it barred designs
% that are fine.
wall = a - h;                          % solid between voids of adjacent cells

% FILLET BOUNDS. Read off the two DiskSelections in buildCrossUnitCell.m
% (:129-148), which decide which vertices each radius is applied to:
%
%   r1 -> disksel1, an annulus w/2 <= rho <= 3w/2 about the origin. That
%         always catches the four REENTRANT corners at (+-w/2, +-w/2), whose
%         distance from the origin is w/sqrt(2).
%   r2 -> disksel2, a thin annulus at rho ~ h/2, which catches the eight
%         ARM-END corners at (+-h/2, +-w/2) and (+-w/2, +-h/2).
%
% An arm end is a rectangle of width w capped by two corners, and both of
% those corners eat into the same end edge, so r2 <= w/2. Along the arm SIDE,
% of length (h-w)/2, a reentrant fillet and an arm-end fillet consume the edge
% from opposite ends, so r1 + r2 <= (h-w)/2. Exceed either and COMSOL's Fillet
% either fails or silently degenerates the corner.
%
% Note what is NOT here: the old code bounded the fillet by the ligament
% (a-h)/2, i.e. by the distance to the cell boundary. That is the same mistake
% as above and is backwards on top of it -- rounding an arm-end corner pulls
% the void AWAY from the boundary, so a fillet can only ever WIDEN the wall,
% never threaten it. Nothing about a fillet depends on a.
rules = { ...
    'inter-cell wall (a-h)',      wall        - minFeature + tol; ...
    'arm width w',                w           - minFeature + tol; ...
    'h > w (is a cross)',         h - w       - tol;              ...
    'arm-end fillet r2 <= w/2',   w/2         - r2 + tol;         ...
    'fillets fit arm side',       (h - w)/2   - (r1 + r2) + tol;  ...
    'positive thickness',         th          - tol};

margin = inf(sz);
why    = repmat({''}, sz);
for k = 1:size(rules, 1)
    mk = rules{k, 2};
    tighter = mk < margin;
    margin(tighter) = mk(tighter);
    why(tighter)    = rules(k, 1);
end

% Non-physical inputs are infeasible by a wide margin, so a graded penalty
% pushes strongly away from them rather than letting a NaN propagate.
bad = ~isfinite(a) | ~isfinite(h) | ~isfinite(w) | ~isfinite(th) | ...
      a <= 0 | h <= 0 | w <= 0;
margin(bad) = -inf;
why(bad)    = {'non-physical dimensions'};

ok      = margin >= 0;
why(ok) = {''};
end
