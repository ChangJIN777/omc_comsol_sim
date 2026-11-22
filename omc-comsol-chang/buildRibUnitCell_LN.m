function [model,P] = buildRibUnitCell_LN(model,P)
%BUILDRIBUNITCELL_LN  'rib' unit cell: spine + centred rib, optional sidewall angle.
%
%   [model,P] = buildRibUnitCell_LN(model,P)
%
%   The single unit-cell builder for celltype 'rib'. The name is historical:
%   this file used to hold the lithium-niobate rib cell, and it now holds the
%   merged builder (the visible-rib cell of comsol_templates/visible_rib_omc.m
%   plus the LN sidewall-angle construction). Its callers are unchanged:
%   runBands, runBands_noSym, runOpticalBand_1D and runOpticalBand_3D. The
%   default material is diamond; set P.beamMat = 'LN' as before for LN.
%
%   GEOMETRY. Spine x in [-a/2, a/2], |y| <= s/2; rib x in [-t/2, t/2],
%   |y| <= L; slab z in [-th/2, th/2]. All widths are at the bottom face.
%     P.d = 0  vertical walls, in-plane fillets r1 (root) and r2 (tip), drawn
%              with the template's feature chain (addRibCell)
%     P.d > 0  walls lean in by th*tan(d) toward the top face, built as
%              tapered prisms (addRibPrisms, the LN construction). No fillets.
%   L follows ribLength: an explicit P.L or P.w wins, otherwise the L(a) rule.
%   The value used is written back to P.L and P.w = 2L.
%
%   DIFFERENCE FROM THE OLD LN CELL. The rib is centred in the cell (the LN
%   cell put half-ribs on x = +-a/2). That is the same lattice shifted by
%   a/2, so the band structure is identical, but the Floquet faces at
%   x = +-a/2 now cut only the spine. The old builder's model parameters a, s,
%   w, t, th, d, airrad, kx, k and its result tables tbl1-tbl4 are no longer
%   created. The runners set a/k/kx/ky themselves, and runOpticalBand_3D
%   creates tbl1-tbl14 itself (it collided with the old tables). Only lbd0,
%   read by the LN refractive index in the optical runners, is still set.
%
%   P FIELDS  see ribDefaultParams / rib_README.md (defaults filled here).
%     P.mbeveny ~= 0 or P.run_optical   only y >= 0 is built; the y = 0 face
%                     goes to 'yboundaries' (the optical runners always treat
%                     y = 0 as a symmetry plane, as the old LN cell did)
%     P.mbevenz ~= 0  z < 0 is cut away (d = 0 only); the z = 0 face(s) are 'ZsymSel'
%     P.run_optical   air cylinder of radius P.airradUC (default 5a) around the
%                     cell, as the LN cell. Optical solvers only: runBands and
%                     runBands_noSym refuse it, because they assign the beam
%                     material with selection.all.
%     Missing P.mbeveny / P.mbevenz are treated as 0.
%
%   PROVIDED TO THE CALLER
%     geom1_xboundaries_bnd   x = +-a/2 faces (beam, plus air when present)
%     geom1_yboundaries_bnd   y = 0 face(s)                  (half-y builds)
%     geom1_ZsymSel           z = 0 face(s)                  (P.mbevenz ~= 0)
%     geom1_beamDomSel, geom1_airDomSel   domain selections  (air: run_optical)
%     geom1_airOuter_bnd      curved outer air surface       (run_optical)
%     P.domSel.beam, P.domSel.air, P.bndSel.Zsym, P.bndSel.airOuter
%     P.zEnd2                 = P.bndSel.airOuter, read by runOpticalBand_3D's
%                             scattering boundary (airTopBoundary)
%     P.xEnd1, P.xEnd2, P.yEnd1, P.zEnd   bndindex lookups
%     P.L, P.w (+ P.ribDerived)           rib length used
%   No flat outer y face exists, so there is no P.yEnd2 and runBands'
%   fixed_bc fallback errors, correctly.

P = ribDefaultParams(P);
validateRibParams(P, 'unitcell');

a  = P.a;   t  = P.t;   s  = P.s;   th = P.th;   d = P.d;
[L, Lsrc] = ribLength(P, a);
P.L = L;
P.w = 2*L;
% Record only what was computed here, so a caller's explicit value stays explicit.
if ~strcmp(Lsrc,'L'), P.ribDerived.L = P.L; end
if ~strcmp(Lsrc,'w'), P.ribDerived.w = P.w; end

mbeveny = 0;  mbevenz = 0;
if isfield(P,'mbeveny'), mbeveny = P.mbeveny; end
if isfield(P,'mbevenz'), mbevenz = P.mbevenz; end
runOpt = isfield(P,'run_optical') && P.run_optical;
symY   = mbeveny ~= 0 || runOpt;
symZ   = mbevenz ~= 0;       % validateRibParams rejects this with d > 0
R      = P.airradUC;
selDelta = 10e-9;

%% Component and geometry. The runners hard-code 'comp1' / 'geom1'.
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label('Rib unit cell');

% Model parameter read by the LN dispersion in runOpticalBand_1D/_3D.
model.param.set('lbd0', sprintf('%.12g[m]', P.lbd0), 'wavelength band edge');

%% Beam solid
if d == 0
    ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
    ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
    addRibCell(ucellWP.geom, 'u', 0, a/2, a/2, t, L, s, P.r1, P.r2, ~symY);
    finTag = 'ext1';
    ext = ucellgeom.feature.create(finTag, 'Extrude');
    ext.set('distance', th);
    ext.selection('input').set({'wp1'});
else
    finTag = addRibPrisms(ucellgeom, 'u', 0, -a/2, a/2, t, L, s, th, d, symY);
end

%% Optional air region (optical runs only), as the LN cell's wp8/rev1
yExt = L;
if runOpt
    wpAir = ucellgeom.feature.create('wpAir', 'WorkPlane');
    wpAir.set('quickplane', 'xz');
    airRect = wpAir.geom.feature.create('airRect', 'Rectangle');
    airRect.set('pos', [-a/2 0]).set('size', [a R]);
    revAir = ucellgeom.feature.create('revAir', 'Revolve');
    revAir.selection('input').set({'wpAir'});
    revAir.set('axis', [-1 0]);
    revAir.set('angle2', 180);          % y >= 0 half (symY is true here)
    cellAir = ucellgeom.feature.create('cellAir', 'Compose');
    cellAir.selection('input').set({finTag 'revAir'});
    cellAir.set('formula', [finTag,' + revAir']);
    finTag = 'cellAir';
    yExt = R;
end

%% z symmetry: remove z < 0 (vertical walls only)
if symZ
    zCut = max(th/2, runOpt*R) + selDelta;
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -zCut);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [-a/2-selDelta, -yExt-selDelta]).set('size', [a+2*selDelta, 2*yExt+2*selDelta]);
    symZExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZExt.set('distance', zCut);
    symZExt.selection('input').set({'symZWP'});
    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set({finTag 'symZPlaneExt'});
    symZComp.set('formula', [finTag,' - symZPlaneExt']);

    % allvertices in this box: the beam's z = 0 face for mechanics, all
    % z = 0 faces (beam + air) when the air is present.
    ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
    ZsymSel.set('xmin', -a/2-selDelta).set('xmax', a/2+selDelta);
    ZsymSel.set('ymin', -yExt-selDelta).set('ymax', yExt+selDelta);
    ZsymSel.set('zmin', -selDelta).set('zmax', selDelta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices');
end

%% Domain selections
beamDomSel = ucellgeom.create('beamDomSel', 'BoxSelection');
beamDomSel.set('xmin', -a/2-selDelta).set('xmax', a/2+selDelta);
beamDomSel.set('ymin', -L-selDelta).set('ymax', L+selDelta);
beamDomSel.set('zmin', -th/2-selDelta).set('zmax', th/2+selDelta);
beamDomSel.set('entitydim', 3).set('condition', 'inside');
if runOpt
    allDomSel = ucellgeom.create('allDomSel', 'BoxSelection');
    allDomSel.set('xmin', -a/2-selDelta).set('xmax', a/2+selDelta);
    allDomSel.set('ymin', -R-selDelta).set('ymax', R+selDelta);
    allDomSel.set('zmin', -R-selDelta).set('zmax', R+selDelta);
    allDomSel.set('entitydim', 3).set('condition', 'inside');
    airDomSel = ucellgeom.create('airDomSel', 'DifferenceSelection');
    airDomSel.set('entitydim', 3).set('add', 'allDomSel').set('subtract', 'beamDomSel');
end

%% Named boundary selections consumed by the runners
% Cumulative selections are created before the boxes that contribute to them.
ucellgeom.selection.create('xboundaries', 'CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
addFaceBox(ucellgeom, 'x_boundary_boxsel_l', 'x', -a/2, selDelta, 'xboundaries');
addFaceBox(ucellgeom, 'x_boundary_boxsel_r', 'x',  a/2, selDelta, 'xboundaries');
if symY
    ucellgeom.selection.create('yboundaries', 'CumulativeSelection');
    ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
    addFaceBox(ucellgeom, 'y_boundary_symmetry', 'y', 0, selDelta, 'yboundaries');
end
if runOpt
    % Curved outer air surface (scattering boundary): small 'intersects'
    % boxes at 45 deg above and below the axis. Only the curved face reaches
    % them; with the z cut only the upper one finds anything.
    ucellgeom.selection.create('airOuter', 'CumulativeSelection');
    ucellgeom.selection('airOuter').label('Cumulative Selection air outer surface');
    for sgn = [1 -1]
        tag = sprintf('airOuterBox%d', (3 - sgn)/2);
        box = ucellgeom.feature.create(tag, 'BoxSelection');
        box.set('entitydim', 2).set('condition', 'intersects');
        box.set('xmin', -selDelta).set('xmax', selDelta);
        box.set('ymin', R/sqrt(2) - 10*selDelta).set('ymax', R/sqrt(2) + 10*selDelta);
        box.set('zmin', sgn*R/sqrt(2) - 10*selDelta).set('zmax', sgn*R/sqrt(2) + 10*selDelta);
        box.set('contributeto', 'airOuter');
    end
end

%% Finalize. run, NOT runAll: runAll leaves the geometry unfinalized and every
% selection / bndindex lookup below would come back empty.
ucellgeom.run;

if symZ
    P.bndSel.Zsym = resolveSel(model, [ucellname,'_ZsymSel'], 2);
end
P.domSel.beam = resolveSel(model, [ucellname,'_beamDomSel'], 3);
if runOpt
    P.domSel.air = resolveSel(model, [ucellname,'_airDomSel'], 3);
    P.bndSel.airOuter = resolveSel(model, [ucellname,'_airOuter_bnd'], 2);
    P.zEnd2 = P.bndSel.airOuter;
end
P.xEnd1 = bndindex(ucellgeom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [ a/2 0 0], [1 0 0]);
if symY
    P.yEnd1 = bndindex(ucellgeom, [0 0 0], [0 1 0]);
end
% bndindex uses p0 only to fix the plane: the bottom face is at z = 0 after the
% z cut, otherwise at z = -th/2.
if symZ
    zProbe = 0;
else
    zProbe = -th/2;
end
P.zEnd = bndindex(ucellgeom, [0 0 zProbe], [0 0 1]);

warnIfEmpty(model, [ucellname,'_xboundaries_bnd'], 2);
if symY,   warnIfEmpty(model, [ucellname,'_yboundaries_bnd'], 2); end
if symZ,   warnIfEmpty(model, [ucellname,'_ZsymSel'], 2); end
if runOpt, warnIfEmpty(model, [ucellname,'_airOuter_bnd'], 2); end
if isempty(P.domSel.beam)
    warning('rib:emptyDomSel', 'geom1_beamDomSel resolved to no domains.');
end
% The runners draw the geometry themselves when P.plotgeom is set.
end

% -------------------------------------------------------------------------

function addFaceBox(geom, tag, axisName, pos, selDelta, contributeTo)
%ADDFACEBOX  'inside' BoxSelection of the faces lying in one axis-normal plane.
box = geom.feature.create(tag, 'BoxSelection');
box.set('entitydim', 2);
box.set([axisName,'min'], pos-selDelta/2).set([axisName,'max'], pos+selDelta/2);
box.set('inputent', 'all');
box.set('condition', 'inside');
box.set('contributeto', contributeTo);
end

% -------------------------------------------------------------------------

function inds = resolveSel(model, tag, dim)
%RESOLVESEL  Entity indices of a geometry selection, as a row vector.
% Two-accessor fallback as in buildCrossStrip's resolveBndSel. Returns []
% rather than throwing. Only meaningful after geom.run.
try
    inds = double(model.selection(tag).entities(dim));
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

function warnIfEmpty(model, tag, dim)
%WARNIFEMPTY  Warn when a selection a solver will hang a BC on is empty.
if isempty(resolveSel(model, tag, dim))
    warning('rib:selectionEmpty', ...
        ['%s resolved to no entities. A periodic, symmetry or scattering ' ...
         'condition using it would be applied to nothing.'], tag);
end
end
