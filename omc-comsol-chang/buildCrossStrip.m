function [model,P] = buildCrossStrip(model,P)
%BUILDCROSSSTRIP  Square-lattice strip (supercell) of N cross unit cells.
%
%   [model,P] = buildCrossStrip(model,P)
%
%   Strip/supercell analogue of buildCrossUnitCell.m, standing in the same
%   relationship to it that buildHoleStrip_3D.m stands to a single hole cell -
%   but on a SQUARE lattice instead of the hexagonal one. The practical
%   difference from buildHoleStrip_3D:
%
%       buildHoleStrip_3D (hexagonal)     buildCrossStrip (square)
%       ------------------------------    ---------------------------------
%       row pitch  sqrt(3)*a/2            row pitch  a          (no sqrt(3))
%       x alternates 0 and +-a/2          x is always 0         (no stagger)
%       half-holes straddle x = +-a/2     nothing is cut by x = +-a/2
%       13 hardcoded circle features      2*P.ncell looped rectangle features
%       extent hardcoded sqrt(3)*4.5*a    extent = P.ncell*a + P.b_wvg
%
%   GEOMETRY
%     A base rectangle with its corner at [-a/2, yLo] and size [a, Ly] is
%     etched with P.ncell cross-shaped voids, all centred on x = 0 and spaced a
%     apart in y, then extruded by P.th. Each void is the union of two crossed
%     rectangles, [h w] and [w h], exactly as in buildCrossUnitCell.
%
%         y_i  = b_wvg + (i - 1/2)*a + b*(i == 1),     i = 1 .. P.ncell
%         yTop = b_wvg + P.ncell*a
%         yLo  = 0     (default)  or  y_1  (P.cutBottomHalfCell ~= 0)
%         Ly   = yTop - yLo
%
%     Ly is the EXTENT (the Rectangle 'size'); yLo and yTop are ABSOLUTE y.
%     They coincide only in the default case, so do not read Ly as the y of the
%     top face.
%
%     With the defaults b = b_wvg = 0 and P.cutBottomHalfCell = 0 this is a
%     plain square lattice: y = 0 and y = yTop land exactly on cell edges, one
%     lattice vector apart, so the two y faces are a legitimate periodic pair.
%
%   HALF-CELL TRUNCATION (P.cutBottomHalfCell ~= 0)
%     The footprint starts at yLo = y_1, the CENTRE of the first cross, so the
%     lower half of the bottom unit cell is gone and the strip is
%     P.ncell - 1/2 cells long. NOTHING MOVES: every cell centre and every void
%     keeps the same absolute y; only the lower extent of the base rectangle,
%     and hence Ly, changes.
%
%     Consequences worth knowing:
%       * The bottom face cuts cell 1 through the middle of its horizontal arm,
%         so it is no longer one rectangle but TWO, each only (a - h)/2 wide in
%         x - the same ligament that already appears back to back across every
%         interior cell boundary. It is still contributed to 'yboundaries' and
%         still returned as P.yEnd1, now as two boundary indices.
%       * The two y faces are then half a lattice vector apart, so they are NO
%         LONGER a periodic pair. Use this only when the y faces carry
%         something other than a Floquet pairing - e.g. runBands with
%         P.fixed_bc = 1 and P.fixed_faces = {'yEnd2'}, which clamps the top
%         face and leaves the cut face free.
%       * The cut creates two new vertices at (+-h/2, yLo), at radius exactly
%         h/2 from the cell-1 centre, i.e. dead centre of the legacy disksel2
%         annulus [h/2 - 5nm, h/2 + 5nm]. So P.r2, a no-op on the arm tips (see
%         KNOWN ISSUE below), DOES round those two corners, which trims each
%         bottom face from (a - h)/2 of flat to (a - h)/2 - r2 of flat plus an
%         r2 arc - the thinnest flat feature in the model, so check the mesh
%         there. checkFilletSelectionBleed reports it under
%         :filletSelectionBleed. P.r2 = 0 avoids it and, wherever
%         :filletR2NoOp is also firing, costs nothing elsewhere.
%
%     h < a is REQUIRED (enforced). At h = a the arms reach the cell boundary
%     and the solid falls into disconnected corner islands.
%
%   TOPOLOGY ASSUMPTION - PLEASE READ
%     This builder uses the half-strip FOOTPRINT of buildHoleStrip_3D (y runs
%     from yLo to yTop, not -Ly/2 to +Ly/2) but deliberately does NOT assume a
%     mirror/symmetry plane at y = 0. buildHoleStrip_3D treats y = 0 as a
%     mirror and therefore emits only ONE y selection; this file emits TWO,
%     at y = yLo and at y = yTop, contributes both to the 'yboundaries'
%     cumulative selection, and returns both P.yEnd1 (y = yLo) and P.yEnd2
%     (y = yTop). Which boundary condition to hang on them - Floquet pair,
%     symmetry, fixed, free - is left entirely to the caller.
%
%     Consequence worth knowing: with b_wvg ~= 0 the lower face is no longer a
%     lattice translation of the y = yTop face, so a Floquet pair across them
%     is no longer the physical square lattice. The faces are still congruent
%     rectangles, so COMSOL will accept the pairing without complaint. With
%     P.cutBottomHalfCell the two faces are not even congruent - see HALF-CELL
%     TRUNCATION above.
%
%   SOLVER TARGET - runBands (MECHANICAL)
%     NO air region is built. buildHoleStrip_3D's air disk (the 'wpair' work
%     plane plus the 90 deg Revolve) is omitted entirely, not merely gated off:
%     runBands_2D pins mbfem.b_domind = 1 and runOpticalBand_* assign materials
%     to hardcoded domain indices, so a second domain would get no material, no
%     physics, and would renumber the diamond. If an optical run is wanted
%     later, add the air region behind an explicit P.add_airDisk flag following
%     buildBoomerangUnitCell.m.
%
%   REQUIRED P FIELDS
%     P.a        lattice constant [m], in BOTH x and y
%     P.h        cross arm length [m]  (must satisfy h < a)
%     P.w        cross arm width  [m]
%     P.th       slab thickness in z [m]
%     P.r1       fillet radius, inner corners (see KNOWN ISSUE below)
%     P.r2       fillet radius, arm tips      (see KNOWN ISSUE below)
%     P.ncell    number of cross cells along the strip (y). P.nperiod is
%                untouched and still means periods along x (= 1 for a strip).
%     P.mbevenz  nonzero -> cut away z < 0 and emit the z symmetry plane
%
%   OPTIONAL P FIELDS
%     P.b        extra y shift applied to cell 1 only, mirroring the role of
%                P.b in buildHoleStrip_3D. Default 0.
%     P.b_wvg    y offset applied to the whole array, widening the gap between
%                y = 0 and the first cell. Default 0.
%     P.cutBottomHalfCell
%                nonzero -> start the footprint at yLo = y_1 instead of y = 0,
%                removing the lower half of the bottom unit cell (see HALF-CELL
%                TRUNCATION above). Default 0, which reproduces the original
%                full-cell footprint exactly.
%     P.plotgeom nonzero -> draw the geometry with mphgeom. Default 0. All
%                plotting is gated: this builder is called once per evaluation
%                inside bayesopt_cross.
%
%   PROVIDED TO THE CALLER
%     named selections   geom1_xboundaries_bnd   (x = -a/2 and x = +a/2)
%                        geom1_yboundaries_bnd   (y = yLo and y = yTop)
%                        geom1_ZsymSel           (z = 0 plane, if P.mbevenz)
%                        geom1_yFixedSel         (y = yTop face ALONE)
%     P.bndSel.Zsym      boundary indices of the z symmetry plane
%     P.bndSel.yFixed    NAME of the y = yTop selection, 'geom1_yFixedSel'.
%                        This is what a fixed BC should be hung on:
%                        fixedBCs.selection.named(P.bndSel.yFixed), the way
%                        runBands.m:373 already does for geom1_ZsymSel. A name
%                        re-resolves on every geometry rebuild; an index list
%                        does not. Distinct from geom1_yboundaries_bnd, which
%                        holds BOTH y faces for the symmetry path.
%     P.bndSel.yFixedInds  the indices that selection resolved to, for
%                        printing/debugging only - do not build a BC from them.
%     P.xEnd1 P.xEnd2    bndindex lookups for the x faces (runBands fixed_bc
%                        FALLBACK path only)
%     P.yEnd1 P.yEnd2    bndindex lookups for the y faces at y = yLo and
%                        y = yTop. P.yEnd1 is TWO boundaries, not one, when
%                        P.cutBottomHalfCell cuts the footprint through a void.
%     P.zEnd             boundaries on the z = 0 plane when P.mbevenz ~= 0,
%                        else on z = -th/2. buildCrossUnitCell probes -th/2
%                        unconditionally, which returns EMPTY once the z < 0
%                        half has been cut away; runOpticalBand_3D then sets a
%                        symmetry plane on nothing. Fixed here.
%     P.Ly               total strip extent in y [m], = P.yHi - P.yLo
%     P.yLo P.yHi        absolute y of the lower and the upper face [m]
%     P.ycentres         1xN cell-centre y coordinates [m]
%     P.ucellname        geometry tag, 'geom1'
%
%   KNOWN ISSUE - LEGACY FILLET SELECTIONS, REPRODUCED DELIBERATELY
%     The two DiskSelections are byte-for-byte the ones in
%     buildCrossUnitCell.m:119-149, translated to each cell centre and given
%     unique tags. Nothing else about them was changed, so that results stay
%     comparable with already-simulated and fabricated designs. They are not
%     what their names suggest:
%
%       disksel1 = annulus [w/2, 1.5*w]                      -> fillet r1
%       disksel2 = annulus [h/2 - 5nm, h/2 + 5nm]            -> fillet r2
%
%     The composed profile has three vertex families, at these radii from the
%     cell centre:
%       inner concave corners (4)  w/sqrt(2)
%       arm-tip corners       (8)  sqrt(h^2 + w^2)/2
%       base-rectangle corners(4)  only at the two ends of the strip
%
%     So disksel2, centred on h/2, misses the tips at sqrt(h^2+w^2)/2 unless
%     w^2/(4h) is under about 5 nm, i.e. roughly w < sqrt(20*h[nm]) nm. For the
%     test_CrossUnitCell.m parameters (a=894, h=863, w=416 nm) the tips sit at
%     479 nm against an annulus at 431.5 +- 5 nm: disksel2 selects NOTHING and
%     r2 is a no-op. Meanwhile disksel1's outer radius 1.5*w = 624 nm exceeds
%     479 nm, so r1 is applied to all 12 corners, not just the 4 inner ones.
%     That is the existing behaviour and it is preserved here on purpose.
%
%     In a strip this also lets one cell's selection reach into its neighbours,
%     whose centres are only a away. checkFilletSelectionBleed below computes
%     that in plain MATLAB before anything is built and WARNS, naming the cells
%     involved. It never changes the geometry.
%
%     Be aware of what that bleed implies once N > 1. The fillets run per cell
%     in sequence, so if cell i's disksel1 encloses cell i+1's tip corners,
%     fil1(i) rounds them BEFORE cell i+1's own disksel1 is evaluated. Each
%     Fillet replaces a vertex with two vertices and an arc, both still near
%     the original radius, so cell i+1's disksel1 can then select those arc
%     endpoints and fil1(i+1) tries to fillet an already-filleted corner. That
%     either double-rounds the corner or makes COMSOL reject the operation
%     ("fillet radius too large"). With the test_CrossUnitCell parameters the
%     bleed warning fires for every cell, so this is not hypothetical: build
%     with P.ncell = 3 and inspect in the GUI before trusting a long strip. The
%     fix would be to narrow the annuli to the analytic radii w/sqrt(2) and
%     sqrt(h^2+w^2)/2, which is deliberately NOT done here because it would
%     change the filleting of existing designs.
%
%   See also BUILDCROSSUNITCELL, BUILDHOLESTRIP_3D, BNDINDEX, RUNBANDS.

%% Read, default and validate the parameters
P = readCrossStripParams(P);

a    = P.a;
h    = P.h;
w    = P.w;
th   = P.th;
r1   = P.r1;
r2   = P.r2;
N    = P.ncell;
b    = P.b;
bwvg = P.b_wvg;
cutHalf = abs(P.cutBottomHalfCell);

% Cell centres and total extent. yTop deliberately does not include b, matching
% buildHoleStrip_3D, where b moves only the first hole. That file never checks
% that the shifted hole still fits inside its footprint; readCrossStripParams
% does, so the two cannot drift apart silently.
ycentres = bwvg + ((1:N) - 0.5)*a;
ycentres(1) = ycentres(1) + b;
yTop = bwvg + N*a;

% Lower extent of the footprint. The truncation moves ONLY where the base
% rectangle starts - the cell centres above are untouched, so every void stays
% at the same absolute y. Ly is from here on a size, never a coordinate: all y
% coordinates below are written as yLo/yTop.
if cutHalf
    yLo = ycentres(1);      % cut through the centre of cell 1
else
    yLo = 0;
end
Ly = yTop - yLo;

P.ycentres = ycentres;
P.Ly       = Ly;
P.yLo      = yLo;
P.yHi      = yTop;

% Pure-MATLAB check of the legacy fillet annuli against the neighbouring cells,
% the strip ends and, when truncating, the two vertices the cut creates. Warn
% only - see KNOWN ISSUE above.
checkFilletSelectionBleed(a, h, w, r1, r2, ycentres, yLo, yTop, cutHalf);

%% Create component
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel  = 'Cross strip';
ucellname   = 'geom1';
P.ucellname = ucellname;

ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% 2D profile: base rectangle minus N cross-shaped voids
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
wpGeom = ucellWP.geom;

baseRect = wpGeom.feature.create('r_base', 'Rectangle');
baseRect.label('Strip footprint');
baseRect.set('base', 'corner');
baseRect.set('pos',  [-a/2 yLo]);
baseRect.set('size', [a Ly]);

% One loop, two rectangles per cell, unique tags. No runCurrent inside the
% loop: buildCrossUnitCell runs the sequence after almost every feature, which
% here would mean O(N) full geometry rebuilds on the bayesopt_cross critical
% path. Everything is created first and built once, below.
subTags = cell(1, 2*N);
for i = 1:N
    tagH = sprintf('c%dh', i);
    tagV = sprintf('c%dv', i);

    recH = wpGeom.feature.create(tagH, 'Rectangle');
    recH.set('base', 'center');
    recH.set('pos',  [0 ycentres(i)]);
    recH.set('size', [h w]);

    recV = wpGeom.feature.create(tagV, 'Rectangle');
    recV.set('base', 'center');
    recV.set('pos',  [0 ycentres(i)]);
    recV.set('size', [w h]);

    subTags{2*i-1} = tagH;
    subTags{2*i}   = tagV;
end

% strjoin on a cellstr, not strcat in a loop: strcat(char,cell) returns a cell,
% and strcat also strips trailing whitespace from char arguments.
composeFormula = ['r_base-', strjoin(subTags, '-')];
composeFeat = wpGeom.feature.create('co1', 'Compose');
composeFeat.set('formula', composeFormula);

ucellgeom.runCurrent;

%% Fillets
% Legacy annulus definitions, reproduced exactly - only the centre offset and
% the feature tags differ from buildCrossUnitCell's addFillet. Per-cell order
% is disksel1, fil1, disksel2, fil2, as in the original.
selection_width = 10e-9;
for i = 1:N
    addCrossFillet(wpGeom, i, ycentres(i), h, w, r1, r2, selection_width);
end
ucellgeom.runCurrent;

%% Extrude
extrudeFeat = ucellgeom.feature.create('ext1', 'Extrude');
extrudeFeat.set('distance', th);
% Input selection set explicitly. The default input silently picks the wrong
% work plane once the sequence holds more than one candidate, which it does as
% soon as the z symmetry block below is added.
extrudeFeat.selection('input').set({'wp1'});
ucellgeom.runCurrent;

%% Symmetry in z
if abs(P.mbevenz)
    symZth = th/2;

    % Footprint derives from Ly, so it tracks P.ncell. buildHoleStrip_3D sizes
    % its ZsymSel from P.airDiskH, which is unrelated to the strip length and
    % only works there by coincidence.
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [-a/2 yLo]).set('size', [a Ly]);
    ucellgeom.runCurrent;

    symZPlaneExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', symZth);
    symZPlaneExt.selection('input').set({'symZWP'});

    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set('ext1');
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', 'ext1 - symZPlaneExt');
    ucellgeom.runCurrent;

    % Beam z-symmetry plane.
    delta = 10e-9;
    ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
    ZsymSel.set('xmin', -a/2-delta).set('xmax', a/2+delta);
    ZsymSel.set('ymin', yLo-delta).set('ymax', yTop+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices');
    ucellgeom.runCurrent;
    % P.bndSel.Zsym is NOT read here. The other builders read it at this point
    % (buildCrossUnitCell.m:104, buildBoomerangUnitCell.m:146, ...), but the
    % selection cannot resolve until the geometry is finalized, so it comes back
    % empty. It is read after the ucellgeom.run below instead.
end

%% Named (cumulative) boundary selections
% These are what the solvers actually consume - runBands.m uses
% 'geom1_xboundaries_bnd' for the Floquet pair and 'geom1_yboundaries_bnd' for
% the y symmetry/antisymmetry pair. Being geometric, they scale with P.ncell
% for free. Contrast buildHoleStrip_3D.m:229-233, whose literal index lists
% (P.xEnd1 = [1 3 4 ...]) were captured from one specific build: they go stale
% whenever the cell count, mbevenz, a fillet radius or the feature list
% changes, and for a builder parameterised on N they are wrong by
% construction. The identical four lists also appear verbatim in
% buildSnowflakeStrip_3D.m, for entirely different geometry.
%
% The cumulative selection nodes are created BEFORE the box selections that
% contribute to them.
sel_delta = 10e-9;

ucellgeom.selection.create('xboundaries','CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');

x_boundary_boxsel_l = ucellgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_boxsel_l.set('entitydim', 2);
x_boundary_boxsel_l.set('xmin', -a/2-sel_delta/2).set('xmax', -a/2+sel_delta/2);
x_boundary_boxsel_l.set('inputent', 'all');
x_boundary_boxsel_l.set('condition', 'inside');
x_boundary_boxsel_l.set('contributeto','xboundaries');

x_boundary_boxsel_r = ucellgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_boxsel_r.set('entitydim', 2);
x_boundary_boxsel_r.set('xmin', a/2-sel_delta/2).set('xmax', a/2+sel_delta/2);
x_boundary_boxsel_r.set('inputent', 'all');
x_boundary_boxsel_r.set('condition', 'inside');
x_boundary_boxsel_r.set('contributeto','xboundaries');

% BOTH y faces, because this strip does not assume a mirror at y = 0. The lower
% one sits at yLo: y = 0 by default, the centre of cell 1 when
% P.cutBottomHalfCell truncates the footprint. In that case it resolves to the
% TWO (a - h)/2-wide ligament faces either side of the halved void, which the
% 'inside' condition picks up just as well as a single full-width face.
ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');

y_boundary_boxsel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
y_boundary_boxsel_b.set('entitydim', 2);
y_boundary_boxsel_b.set('ymin', yLo-sel_delta/2).set('ymax', yLo+sel_delta/2);
y_boundary_boxsel_b.set('inputent', 'all');
y_boundary_boxsel_b.set('condition', 'inside');
y_boundary_boxsel_b.set('contributeto','yboundaries');

y_boundary_boxsel_t = ucellgeom.feature.create('y_boundary_boxsel_t', 'BoxSelection');
y_boundary_boxsel_t.set('entitydim', 2);
y_boundary_boxsel_t.set('ymin', yTop-sel_delta/2).set('ymax', yTop+sel_delta/2);
y_boundary_boxsel_t.set('inputent', 'all');
y_boundary_boxsel_t.set('condition', 'inside');
y_boundary_boxsel_t.set('contributeto','yboundaries');

%% Named selection for the FIXED boundary condition: the y = yTop face alone
% Separate node from y_boundary_boxsel_t on purpose. That one contributes to
% 'yboundaries', which deliberately holds BOTH y faces because the symmetry /
% antisymmetry path in runBands.m:362,:367 needs the pair; a fixed condition
% must hit the top face ONLY, so it needs a selection of its own.
%
% This is the ZsymSel pattern (a standalone geometry BoxSelection exposed as
% geom1_<tag> and consumed with selection.named), which is the one mechanism in
% this file already proven to reach the physics intact: runBands.m:373 has been
% hanging the z symmetry condition on 'geom1_ZsymSel' all along. Unlike a
% captured index list it re-resolves on every geometry rebuild, so it cannot go
% stale when P.ncell, P.mbevenz or a fillet radius changes.
%
% The z span covers BOTH cases in one box: z in [0, th/2] when P.mbevenz cut
% the lower half away, z in [-th/2, th/2] when it did not.
yFixedSel = ucellgeom.create('yFixedSel', 'BoxSelection');
yFixedSel.set('xmin', -a/2-sel_delta).set('xmax', a/2+sel_delta);
yFixedSel.set('ymin', yTop-sel_delta).set('ymax', yTop+sel_delta);
yFixedSel.set('zmin', -th/2-sel_delta).set('zmax', th/2+sel_delta);
yFixedSel.set('entitydim', 2).set('condition', 'allvertices');

% Run the whole sequence so the named selections resolve and the FINALIZED
% geometry exists before this function returns.
%
% THIS MUST BE run, NOT runAll. runAll builds the geometry objects of the
% sequence but does not finalize (it does not run the 'fin' feature), and every
% query below reads the finalized geometry: measured on COMSOL 6.3, after
% runAll the sequence reports getNVertices = 0, getNBoundaries = 0,
% getNDomains = 0, so bndindex returns [] for all four faces and the
% cumulative selections resolve to nothing. After run it reports 224 vertices,
% 114 boundaries and 1 domain and every lookup below resolves. The repo-wide
% runAll idiom (buildCrossUnitCell.m:71, buildBoomerangUnitCell.m:264,
% buildHoleStrip_3D.m:180) survives only because the callers later reach
% runBands.m:203 / runBands_2D, which call geom('geom1').run - too late for any
% P.xEnd*/P.yEnd*/P.zEnd the builder itself returns. buildHoleStrip_3D.m:184
% has that call sitting commented out for the same reason.
ucellgeom.run;

%% Selections that need the FINALIZED geometry
% Z symmetry plane, deferred from the mbevenz block above for the reason given
% there. Same accessor pair as checkNamedBndSel, for the same version-to-version
% reason.
if abs(P.mbevenz)
    P.bndSel.Zsym = resolveBndSel(model, [ucellname,'_ZsymSel']);
end

% The fixed-BC selection. The caller gets the NAME, not the indices, so that
% runBands can do fixedBCs.selection.named(P.bndSel.yFixed) and never has to
% know anything about this geometry. P.bndSel.yFixedInds is informational (the
% geomOnly print) - do not build a boundary condition out of it.
P.bndSel.yFixed     = [ucellname,'_yFixedSel'];
P.bndSel.yFixedInds = resolveBndSel(model, P.bndSel.yFixed);
if isempty(P.bndSel.yFixedInds)
    warning('buildCrossStrip:yFixedSelEmpty', ...
        ['%s resolved to no boundaries, so a fixed condition using it would ' ...
         'constrain nothing. The box spans y = [%g, %g] nm; the y = yTop face ' ...
         'should be at %g nm.'], P.bndSel.yFixed, ...
        (yTop-sel_delta)*1e9, (yTop+sel_delta)*1e9, yTop*1e9);
end

%% Index-based selections (runBands fixed_bc path, and P.zEnd)
P.xEnd1 = bndindex(ucellgeom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [ a/2 0 0], [1 0 0]);

P.yEnd1 = bndindex(ucellgeom, [0 yLo  0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [0 yTop 0], [0 1 0]);

% bndindex uses p0 only to define a plane (bndindex.m:39-45), so only the z
% coordinate matters here. The z < 0 half is gone when P.mbevenz is nonzero, so
% the surviving bottom face is at z = 0, not at z = -th/2.
if abs(P.mbevenz)
    zProbe = 0;
else
    zProbe = -th/2;
end
P.zEnd = bndindex(ucellgeom, [0 0 zProbe], [0 0 1]);

% Warn-only cross-check of the named selections against the index lookups.
checkNamedBndSel(model, ucellname, 'xboundaries', [P.xEnd1 P.xEnd2]);
checkNamedBndSel(model, ucellname, 'yboundaries', [P.yEnd1 P.yEnd2]);

%% Optional geometry plot
if P.plotgeom
    figure;
    mphgeom(model);
end
end

% -------------------------------------------------------------------------

function P = readCrossStripParams(P)
%READCROSSSTRIPPARAMS Fill in defaults and validate, before COMSOL sees anything.
%
% Catching these here matters: a bad size or a zero extrude distance surfaces
% from COMSOL as an opaque Java geometry error several calls deeper.

required = {'a','h','w','th','r1','r2','ncell','mbevenz'};
missing  = required(~isfield(P, required));
if ~isempty(missing)
    error('buildCrossStrip:missingFields', ...
        'Missing required P field(s): %s.', strjoin(missing, ', '));
end

validateattributes(P.a,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.a');
validateattributes(P.h,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.h');
validateattributes(P.w,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.w');
validateattributes(P.th, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.th');
validateattributes(P.r1, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r1');
validateattributes(P.r2, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r2');
validateattributes(P.ncell, {'numeric'}, {'scalar','real','finite','integer','positive'}, mfilename, 'P.ncell');
validateattributes(P.mbevenz, {'numeric','logical'}, {'scalar','real','finite'}, mfilename, 'P.mbevenz');

% Hard error: at h >= a the arms reach the cell boundary and the remaining
% solid breaks into disconnected corner islands. COMSOL builds it happily.
if P.h >= P.a
    error('buildCrossStrip:armTooLong', ...
        ['P.h = %g nm must be strictly less than P.a = %g nm. At h >= a the ' ...
         'cross arms reach the cell boundary and the solid separates into ' ...
         'disconnected corner islands.'], P.h*1e9, P.a*1e9);
end

if P.w >= P.h
    warning('buildCrossStrip:armWidthExceedsLength', ...
        ['P.w = %g nm is not smaller than P.h = %g nm, so the cross void is ' ...
         'degenerate (a square at w = h). The fillet selections assume a ' ...
         'genuine plus shape.'], P.w*1e9, P.h*1e9);
end

if ~isfield(P,'b'),        P.b = 0;        end
if ~isfield(P,'b_wvg'),    P.b_wvg = 0;    end
if ~isfield(P,'plotgeom'), P.plotgeom = 0; end

% Default 0 = the original full-cell footprint, so existing callers
% (bayesopt_cross, cross_optimize_sweep_diamond, ...) are unaffected.
if ~isfield(P,'cutBottomHalfCell'), P.cutBottomHalfCell = 0; end

validateattributes(P.b,     {'numeric'}, {'scalar','real','finite'}, mfilename, 'P.b');
validateattributes(P.b_wvg, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.b_wvg');
validateattributes(P.cutBottomHalfCell, {'numeric','logical'}, ...
    {'scalar','real','finite','binary'}, mfilename, 'P.cutBottomHalfCell');

% Cell 1 carries the extra shift b while Ly does not, so it is the only cell
% that can be pushed outside the footprint. buildHoleStrip_3D has exactly this
% desynchronisation and never checks for it.
y1 = P.b_wvg + P.a/2 + P.b;
if y1 - P.h/2 < 0
    error('buildCrossStrip:firstCellOutsideFootprint', ...
        ['P.b = %g nm pushes the first cross void below y = 0 (its lower ' ...
         'edge lands at %g nm). Increase P.b_wvg or reduce |P.b|.'], ...
        P.b*1e9, (y1 - P.h/2)*1e9);
end
if y1 + P.h/2 > P.b_wvg + P.ncell*P.a
    error('buildCrossStrip:firstCellOutsideFootprint', ...
        'P.b = %g nm pushes the first cross void past the end of the strip.', ...
        P.b*1e9);
end

% Half-cell truncation. The cut lands on y1 by construction, so it always
% bisects cell 1; what it must not do is leave nothing behind or reach into
% cell 2. Only P.b can take it there (b moves cell 1 but not the strip ends),
% which is the same desynchronisation the two checks above guard.
if P.cutBottomHalfCell
    yTop = P.b_wvg + P.ncell*P.a;
    if y1 >= yTop
        error('buildCrossStrip:truncationLeavesNothing', ...
            ['P.cutBottomHalfCell truncates the footprint at the centre of ' ...
             'cell 1 (y = %g nm), which is at or past the end of the strip ' ...
             '(y = %g nm), so no solid would remain.'], y1*1e9, yTop*1e9);
    end
    if P.ncell > 1
        y2lower = P.b_wvg + 1.5*P.a - P.h/2;    % lower edge of cell 2's void
        if y1 >= y2lower
            error('buildCrossStrip:truncationCutsSecondCell', ...
                ['P.cutBottomHalfCell truncates the footprint at y = %g nm, ' ...
                 'which is already inside the SECOND cross void (its lower ' ...
                 'edge is at %g nm). The truncation is only meant to halve ' ...
                 'cell 1; reduce P.b.'], y1*1e9, y2lower*1e9);
        end
    end
end
end

% -------------------------------------------------------------------------

function checkFilletSelectionBleed(a, h, w, r1, r2, ycentres, yLo, yTop, cutHalf)
%CHECKFILLETSELECTIONBLEED Warn when a legacy fillet annulus reaches another cell.
%
% The two DiskSelections are reproduced from buildCrossUnitCell unchanged, and
% the outer one is broad: [w/2, 1.5*w]. In a single cell the only other
% vertices are the square corners at a/sqrt(2); in a strip the neighbouring
% cells sit only a away, so an annulus can enclose their vertices and the
% fillet then rounds corners belonging to a different cell.
%
% This is a pure-MATLAB pre-build check. It warns and returns; it never alters
% the geometry, because the whole point of the legacy radii is bit-for-bit
% comparability with earlier simulations.

N = numel(ycentres);

% The 12 vertices of the composed profile around one cell centre: 8 arm-tip
% corners and 4 inner concave corners.
vx = [ h/2  h/2 -h/2 -h/2  w/2  w/2 -w/2 -w/2  w/2  w/2 -w/2 -w/2].';
vy = [ w/2 -w/2  w/2 -w/2  h/2 -h/2  h/2 -h/2  w/2 -w/2  w/2 -w/2].';

allX = repmat(vx, N, 1);
allY = reshape(vy + ycentres(:).', [], 1);
owner = repmat(1:N, numel(vx), 1);
owner = owner(:);

% Plus the four outer corners of the strip footprint, owner 0.
allX  = [allX;  -a/2;  a/2; -a/2;  a/2];
allY  = [allY;   yLo;  yLo; yTop; yTop];
owner = [owner;    0;    0;    0;    0];

% And, when the bottom half of cell 1 is cut away, the two vertices the cut
% leaves where it crosses the horizontal arm of that void: (+-h/2, yLo). They
% sit at radius exactly h/2 from the cell-1 centre, which is the middle of the
% legacy disksel2 annulus, so they are the one case where r2 stops being a
% no-op. Owner -1 so they read as foreign to every cell and get reported.
if cutHalf
    allX  = [allX;  -h/2;  h/2];
    allY  = [allY;   yLo;  yLo];
    owner = [owner;   -1;   -1];
end

% The legacy annuli, in the same order the fillets are applied.
sel_width = 10e-9;
annuli = {};
if r1 > 0
    annuli{end+1} = struct('name','disksel1 (r1)', 'rin', w/2,               'rout', w*(3/2));
end
if r2 > 0
    annuli{end+1} = struct('name','disksel2 (r2)', 'rin', h/2-sel_width/2,   'rout', h/2+sel_width/2);
end

msgs = {};
for k = 1:numel(annuli)
    A = annuli{k};
    for i = 1:N
        foreign = owner ~= i;
        d = hypot(allX(foreign), allY(foreign) - ycentres(i));
        hitOwner = owner(foreign);
        inAnnulus = (d >= A.rin) & (d <= A.rout);
        if any(inAnnulus)
            culprits = unique(hitOwner(inAnnulus));
            cellHits = culprits(culprits > 0);
            edgeHit  = any(culprits == 0);
            cutHit   = any(culprits == -1);
            parts = {};
            if ~isempty(cellHits)
                parts{end+1} = ['cell(s) ', num2str(cellHits(:).')]; %#ok<AGROW>
            end
            if edgeHit
                parts{end+1} = 'the strip end corners'; %#ok<AGROW>
            end
            if cutHit
                parts{end+1} = ['the two vertices left by the ' ...
                    'P.cutBottomHalfCell cut at y = ', ...
                    num2str(yLo*1e9,'%.1f'), ' nm']; %#ok<AGROW>
            end
            msgs{end+1} = sprintf('  %s on cell %d (annulus [%.1f, %.1f] nm) reaches %s', ...
                A.name, i, A.rin*1e9, A.rout*1e9, strjoin(parts, ' and ')); %#ok<AGROW>
        end
    end
end

if ~isempty(msgs)
    warning('buildCrossStrip:filletSelectionBleed', ...
        ['The legacy fillet DiskSelections reach beyond their own cell, so a ' ...
         'fillet will round corners that belong elsewhere:\n%s\n' ...
         'These annuli are reproduced from buildCrossUnitCell on purpose (see ' ...
         'the KNOWN ISSUE block in buildCrossStrip.m) and have NOT been ' ...
         'adjusted. Reduce P.w relative to P.a, or accept the extra fillets.'], ...
        strjoin(msgs, newline));
end

% Separately: flag the case where the intended tip fillet catches nothing,
% which is the single-cell behaviour and is easy to miss.
if r2 > 0
    Rtip = sqrt(h^2 + w^2)/2;
    if Rtip < h/2 - sel_width/2 || Rtip > h/2 + sel_width/2
        warning('buildCrossStrip:filletR2NoOp', ...
            ['disksel2 spans [%.1f, %.1f] nm but the arm-tip corners sit at ' ...
             '%.1f nm, so it selects nothing and P.r2 = %g nm has no effect. ' ...
             'This reproduces buildCrossUnitCell exactly; it is reported, not ' ...
             'corrected.'], ...
            (h/2-sel_width/2)*1e9, (h/2+sel_width/2)*1e9, Rtip*1e9, r2*1e9);
    end
end
end

% -------------------------------------------------------------------------

function addCrossFillet(wpGeom, idx, yc, h, w, r1, r2, selection_width)
%ADDCROSSFILLET Legacy cross fillets, translated to one cell centre.
%
% Byte-for-byte the selections from buildCrossUnitCell.m:119-149. The ONLY
% differences are structurally required by the strip: posy is the cell centre
% instead of 0, and every feature tag carries the cell index. Radii, the rin/r
% pairs, the 'allvertices' condition and the r1/r2 assignment are unchanged.
%
% A zero radius is skipped: COMSOL's Fillet errors on radius 0.

disksel1_label = sprintf('c%d_disksel1', idx);
disksel2_label = sprintf('c%d_disksel2', idx);
fillet1_label  = sprintf('c%d_fil1',     idx);
fillet2_label  = sprintf('c%d_fil2',     idx);

if r1 > 0
    disksel1 = wpGeom.feature.create(disksel1_label, 'DiskSelection');
    disksel1.set('entitydim', 0);
    disksel1.set('posx', 0);
    disksel1.set('posy', yc);
    disksel1.set('r', w*(3/2));
    disksel1.set('rin', w/2);
    disksel1.set('condition', 'allvertices');

    fil1 = wpGeom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
end

if r2 > 0
    disksel2 = wpGeom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', yc);
    disksel2.set('posx', 0);
    disksel2.set('r', h/2+selection_width/2);
    disksel2.set('rin', h/2-selection_width/2);
    disksel2.set('condition', 'allvertices');

    fil2 = wpGeom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end
end

% -------------------------------------------------------------------------

function inds = resolveBndSel(model, tag)
%RESOLVEBNDSEL Boundary indices of a geometry selection, as a row vector.
%
% Two accessors because they differ across COMSOL versions, the same pair (and
% the same order) checkNamedBndSel uses. Returns [] rather than throwing if
% neither works: the caller decides whether an empty selection is fatal.
%
% Only meaningful AFTER the geometry has been finalized with geom.run - before
% that every selection resolves to nothing. See the comment on the run call.

try
    inds = double(model.selection(tag).entities(2));
catch
    try
        inds = double(model.selection(tag).inputEntities());
    catch
        inds = [];
    end
end
inds = inds(:)';
end

% -------------------------------------------------------------------------

function checkNamedBndSel(model, ucellname, selname, expected)
%CHECKNAMEDBNDSEL Warn if a cumulative selection misses the expected faces.
%
% Ported from buildBoomerangUnitCell.m. selname is a geometry cumulative
% selection, e.g. 'xboundaries'; expected is the boundary indices the bndindex
% lookups found for that pair.
%
% Reading the resolved entities is wrapped in try/catch: the accessor differs
% across COMSOL versions, and a check that cannot run must not take the
% geometry build down with it.

tag = [ucellname,'_',selname,'_bnd'];
try
    got = double(model.selection(tag).entities(2));
catch
    try
        got = double(model.selection(tag).inputEntities());
    catch err
        warning('buildCrossStrip:selectionUnreadable', ...
            ['Could not read selection %s (%s), so it was not verified ' ...
             'against the bndindex lookup.'], tag, err.message);
        return
    end
end

got = sort(got(:)).';
expected = sort(double(expected(:))).';

if isempty(got)
    warning('buildCrossStrip:selectionEmpty', ...
        ['%s resolved to no boundaries. A periodic or symmetry condition ' ...
         'using it would be applied to nothing.'], tag);
elseif ~isequal(got, expected)
    warning('buildCrossStrip:selectionMismatch', ...
        ['%s resolved to boundaries [%s] but bndindex gives [%s]. Check the ' ...
         'geometry before trusting the boundary conditions.'], ...
        tag, num2str(got), num2str(expected));
end
end
