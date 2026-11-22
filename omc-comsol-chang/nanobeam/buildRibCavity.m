function [model,P] = buildRibCavity(model,P)
%BUILDRIBCAVITY  Rib nanobeam cavity geometry for the FEM pipeline.
%
%   [model,P] = buildRibCavity(model,P)
%
%   RunNanobeamFEM dispatches here for P.celltype = 'rib', after
%   createRibCavityGeom has filled P.ribGeomHalf. It mirrors
%   BuildNanobeamFEMGen.m: same component / geometry / feature tags ('comp',
%   'beam', beamWP, beamHoles, cyl_cut_wp, air_cyl, beamHolesAir, symZWP, ...)
%   and the same P.domSel / P.bndSel outputs, so SetupNanobeamFEM,
%   SolveNanobeamFEM and CalcGOM consume it unchanged.
%
%   GEOMETRY. The x >= 0, y >= 0 quarter of the symmetric cavity of
%   comsol_templates/visible_rib_omc.m, plus z >= 0 when z symmetry is on (same
%   switch as BuildNanobeamFEMGen). The template kept the mirror-image octant
%   x <= 0, y <= 0, which is equivalent. Each cell is drawn by addRibCell
%   with its spine segment [xL, xR] from P.ribGeomHalf. The x < 0 half of the
%   centre cell is clipped with a mask. The beam is centred on z = 0
%   (work plane at -th/2).
%   SIDEWALL ANGLE. With P.d > 0 (radians) each cell is built instead from the
%   tapered prisms of addRibPrisms (walls lean in by th*tan(d) toward the top
%   face, no fillets) and the union is clipped to x >= 0 by a block. z symmetry
%   is then rejected, because the slab is not mirror-symmetric in z.
%
%   AIR CYLINDER. Radius P.airrad, length beamLenHalf (the template's cylinder
%   spans the beam exactly), axis at z = P.airZOffset. The default 60 nm
%   reproduces the template, whose cylinder axis sits 60 nm above its effective
%   symmetry plane; see ribDefaultParams, "WHY 110 nm". Set 0 for a
%   centred cylinder. Built only when P.solveOpt, as in BuildNanobeamFEMGen,
%   whose Revolve set-up is copied verbatim apart from the axis offset.
%
%   REQUIRED P FIELDS
%     P.ribGeomHalf, P.ribGeomCols, P.beamLenHalf   createRibCavityGeom
%     P.solveMech, P.solveOpt, P.mevenz, P.oevenz    as BuildNanobeamFEMGen
%   PROVIDED TO THE CALLER
%     P.geomname = 'beam', P.xc = 0
%     P.domSel.beam, P.domSel.cyl (if air)
%     P.bndSel.beamXsym, beamXend, beamYsym, beamZsym ([] without z cut)
%     P.bndSel.cylXsym, cylXend, cylYsym, cylZsym, cylCurv (if air)
%     named selection 'beam_beamSel'. CalcGOM builds 'beamBndsAll' from the
%     geometry feature 'beamSel', so that tag must not change.
%   The mechanical end condition (P.mEndBC) is applied by SetupNanobeamFEM on
%   P.bndSel.beamXend.
%
%   NOT BUILT: mechanical PML (P.solveMechPML), asymmetric cavity (P.asymCav).

if exist('addRibCell','file') ~= 2 || exist('addRibPrisms','file') ~= 2
    error('rib:pathMissing', ...
        ['addRibCell.m / addRibPrisms.m (omc-comsol-chang root) are not on the path. ' ...
         'Add it with addpath(<root>,''-end''), as test_nanobeamRibFEM.m does.']);
end
P = ribDefaultParams(P);
validateRibParams(P, 'cavity');
required = {'ribGeomHalf','ribGeomCols','beamLenHalf','solveMech','solveOpt','mevenz','oevenz'};
missing = required(~isfield(P, required));
if ~isempty(missing)
    error('rib:missingFields', ...
        'Missing P field(s): %s. Run createRibCavityGeom first.', strjoin(missing, ', '));
end
if isfield(P,'solveMechPML') && P.solveMechPML
    error('rib:noPML', 'buildRibCavity does not build a mechanical PML.');
end
if isfield(P,'asymCav') && P.asymCav
    error('rib:noAsym', 'buildRibCavity builds symmetric cavities only.');
end

C       = P.ribGeomCols;
G       = P.ribGeomHalf;
ncell   = size(G,1);
thi     = P.th;
beamLen = P.beamLenHalf;
LMax    = max(G(:,C.L));
P.xc    = 0;

%% Component and geometry (BuildNanobeamFEMGen tags)
comp = model.modelNode.create('comp');
comp.label('Nanobeam FEM simulation');
P.geomname = 'beam';
beamgeom = model.geom.create(P.geomname, 3);
beamgeom.label('Rib nanobeam geometry');

%% Cells, y >= 0 only (the y = 0 plane is a symmetry plane)
selDelta = 10e-9;
finBeamTag = 'beamHoles';
if P.d == 0
    % Vertical walls: template feature chain in a work plane, then extrude.
    beamWP = beamgeom.feature.create('beamWP', 'WorkPlane');
    beamWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -thi/2);
    wpGeom = beamWP.geom;
    cellTags = cell(1, ncell);
    for k = 1:ncell
        xck = G(k,C.x);
        cellTags{k} = addRibCell(wpGeom, sprintf('c%02d',k), xck, ...
            xck - G(k,C.xL), G(k,C.xR) - xck, G(k,C.t), G(k,C.L), G(k,C.s), ...
            P.r1, P.r2, false);
    end

    % Clip to x >= 0. Only the centre cell (P.ribatctr = 1) crosses x = 0.
    halfMask = wpGeom.feature.create('halfMask', 'Rectangle');
    halfMask.set('pos', [0 -selDelta]).set('size', [beamLen, LMax + 2*selDelta]);

    beamComp = wpGeom.feature.create('beamComp', 'Compose');
    beamComp.selection('input').set([cellTags, {'halfMask'}]);
    beamComp.set('formula', ['(', strjoin(cellTags, ' + '), ') * halfMask']);
    beamComp.set('intbnd', false);

    beamHoles = beamgeom.feature.create(finBeamTag, 'Extrude');
    beamHoles.set('distance', thi);
    beamHoles.selection('input').set({'beamWP'});
else
    % Sidewall angle: one tapered spine + rib prism set per cell
    % (addRibPrisms, the LN construction), united, then clipped to x >= 0 by a
    % block. The result keeps the tag 'beamHoles' so everything below is shared.
    cellTags = cell(1, ncell);
    for k = 1:ncell
        cellTags{k} = addRibPrisms(beamgeom, sprintf('c%02d',k), G(k,C.x), ...
            G(k,C.xL), G(k,C.xR), G(k,C.t), G(k,C.L), G(k,C.s), thi, P.d, true);
    end
    halfBox = beamgeom.feature.create('halfBox', 'Block');
    halfBox.set('base', 'corner');
    halfBox.set('pos', [0, -selDelta, -thi/2-selDelta]);
    halfBox.set('size', [beamLen, LMax + 2*selDelta, thi + 2*selDelta]);
    beamHoles = beamgeom.feature.create(finBeamTag, 'Compose');
    beamHoles.selection('input').set([cellTags, {'halfBox'}]);
    beamHoles.set('formula', ['(', strjoin(cellTags, ' + '), ') * halfBox']);
    beamHoles.set('intbnd', false);
end

totLen = beamLen;
maxWid = 2*LMax;
maxThi = thi;

%% Air cylinder (BuildNanobeamFEMGen.m:252-276, plus the z offset of the axis)
hasAir = isfield(P,'airrad') && P.solveOpt;
R    = P.airrad;
zAir = P.airZOffset;
if hasAir
    cyl_cut_wp = beamgeom.feature.create('cyl_cut_wp', 'WorkPlane');
    cyl_cut_wp.set('planetype', 'quick').set('quickplane', 'xz');
    cyl_cut = cyl_cut_wp.geom.feature.create('cyl_cut', 'Rectangle');
    cyl_cut.set('type','solid').set('pos',[0,zAir]).set('size',[totLen,R]);

    air_cyl = beamgeom.feature.create('air_cyl', 'Revolve');
    air_cyl.selection('input').set('cyl_cut_wp');
    air_cyl.set('angle1','-180');
    air_cyl.set('axis',[1,0]).set('pos',[0,zAir]).set('angle2','0');

    beamHolesAir = beamgeom.feature.create('beamHolesAir', 'Compose');
    beamHolesAir.selection('input').set({finBeamTag 'air_cyl'});
    beamHolesAir.set('formula', [finBeamTag,' + air_cyl']);
    finBeamTag = 'beamHolesAir';

    maxWid = max(maxWid, 2*R);
    maxThi = max(maxThi, 2*(R + abs(zAir)));
end

%% z symmetry, same switch as BuildNanobeamFEMGen.m:330-333
symZOn = strcmp(P.xsect,'rect') && ...
        ((P.solveMech && ~P.solveOpt && abs(P.mevenz)) || ...
         (P.solveOpt && ~P.solveMech && abs(P.oevenz)) || ...
         (P.solveMech && P.solveOpt && abs(P.mevenz) && abs(P.oevenz)));
if symZOn && P.d > 0
    error('rib:zSymWithSidewall', ...
        ['z symmetry is on (mevenz/oevenz with solveMech/solveOpt), but with ' ...
         'P.d ~= 0 the slab is not mirror-symmetric in z. Set P.mevenz = 0 ' ...
         'and/or P.oevenz = 0.']);
end
if symZOn
    symZWP = beamgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -maxThi/2);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [0 0]).set('size', [totLen maxWid/2]);
    symZPlaneExt = beamgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', maxThi/2);
    symZPlaneExt.selection('input').set({'symZWP'});
    symZComp = beamgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set({finBeamTag 'symZPlaneExt'});
    symZComp.set('formula', [finBeamTag,' - symZPlaneExt']);
end
zLoBeam = -thi/2*(~symZOn);

%% Domain selections
d = 10e-9;
beamSel = beamgeom.create('beamSel', 'BoxSelection');
beamSel.set('xmin', -d).set('xmax', beamLen+d);
beamSel.set('ymin', -d).set('ymax', LMax+d);
beamSel.set('zmin', -thi/2-d).set('zmax', thi/2+d);
beamSel.set('entitydim', 3).set('condition', 'inside');
if hasAir
    beamCylSel = beamgeom.create('beamCylSel', 'BoxSelection');
    beamCylSel.set('xmin', -d).set('xmax', beamLen+d);
    beamCylSel.set('ymin', -d).set('ymax', R+d);
    beamCylSel.set('zmin', zAir-R-d).set('zmax', zAir+R+d);
    beamCylSel.set('entitydim', 3).set('condition', 'inside');
    cylSel = beamgeom.create('cylSel', 'DifferenceSelection');
    cylSel.set('entitydim', 3).set('add', 'beamCylSel').set('subtract', 'beamSel');
end

%% Boundary selections (flat: allvertices; curved: intersects)
d = 5e-9;
addBox(beamgeom, 'beamXsymSel', [-d d], [-d LMax+d], [zLoBeam-d thi/2+d], 'allvertices');
addBox(beamgeom, 'beamXendSel', [beamLen-d beamLen+d], [-d LMax+d], [zLoBeam-d thi/2+d], 'allvertices');
addBox(beamgeom, 'beamYsymSel', [-d beamLen+d], [-d d], [zLoBeam-d thi/2+d], 'allvertices');
if symZOn
    % The whole beam footprint at z = 0. 'allvertices' excludes the air's z = 0
    % face, which reaches y = airrad. BuildNanobeamFEMGen's tiny 'intersects'
    % box at y = min(w)/2 assumes a solid slab and would land in air here.
    addBox(beamgeom, 'beamZsymSel', [-d beamLen+d], [-d LMax+d], [-d d], 'allvertices');
end
if hasAir
    zLoCyl = (zAir - R)*(~symZOn);
    zHiCyl = zAir + R;
    addBox(beamgeom, 'cylXsymSel', [-d d], [-d R+d], [zLoCyl-d zHiCyl+d], 'allvertices');
    addBox(beamgeom, 'cylXendSel', [beamLen-d beamLen+d], [-d R+d], [zLoCyl-d zHiCyl+d], 'allvertices');
    addBox(beamgeom, 'cylYsymSel', [-d beamLen+d], [-d d], [zLoCyl-d zHiCyl+d], 'allvertices');
    if symZOn
        addBox(beamgeom, 'cylZsymSel', [-d beamLen+d], [-d R+d], [-d d], 'allvertices');
    end
    addBox(beamgeom, 'cylCurvSel', [d 2*d], R/sqrt(2)+[-10*d 10*d], ...
        zAir + R/sqrt(2)+[-10*d 10*d], 'intersects');
end

%% Finalize, then resolve. Selections read before geom.run come back empty.
beamgeom.run;
gn = P.geomname;
P.domSel.beam = resolveSel(model, [gn,'_beamSel'], 3);
if hasAir
    P.domSel.cyl = resolveSel(model, [gn,'_cylSel'], 3);
end
P.bndSel.beamXsym = resolveSel(model, [gn,'_beamXsymSel'], 2);
P.bndSel.beamXend = resolveSel(model, [gn,'_beamXendSel'], 2);
P.bndSel.beamYsym = resolveSel(model, [gn,'_beamYsymSel'], 2);
P.bndSel.beamZsym = [];
if symZOn
    P.bndSel.beamZsym = resolveSel(model, [gn,'_beamZsymSel'], 2);
end
if hasAir
    P.bndSel.cylXsym = resolveSel(model, [gn,'_cylXsymSel'], 2);
    P.bndSel.cylXend = resolveSel(model, [gn,'_cylXendSel'], 2);
    P.bndSel.cylYsym = resolveSel(model, [gn,'_cylYsymSel'], 2);
    P.bndSel.cylZsym = [];
    if symZOn
        P.bndSel.cylZsym = resolveSel(model, [gn,'_cylZsymSel'], 2);
    end
    P.bndSel.cylCurv = resolveSel(model, [gn,'_cylCurvSel'], 2);
end

% Every face selection the setup indexes into must be non-empty, or the
% matching BC is silently applied to nothing.
checkNames = {'beamXsym','beamXend','beamYsym'};
if symZOn, checkNames{end+1} = 'beamZsym'; end
if hasAir, checkNames = [checkNames, {'cylXsym','cylXend','cylYsym','cylCurv'}]; end
for k = 1:numel(checkNames)
    if isempty(P.bndSel.(checkNames{k}))
        warning('rib:emptyBndSel', ...
            'P.bndSel.%s resolved to no boundaries.', checkNames{k});
    end
end
if isempty(P.domSel.beam)
    warning('rib:emptyDomSel', 'P.domSel.beam resolved to no domains.');
end

airStr = '';
if hasAir
    airStr = ', air cylinder';
end
disp(['Geometry created - rib nanobeam, ', num2str(ncell), ...
      ' cells in half beam', airStr]);
end

% -------------------------------------------------------------------------

function addBox(geom, tag, xr, yr, zr, condition)
%ADDBOX  Boundary (entitydim 2) BoxSelection over [xr] x [yr] x [zr].
box = geom.create(tag, 'BoxSelection');
box.set('xmin', xr(1)).set('xmax', xr(2));
box.set('ymin', yr(1)).set('ymax', yr(2));
box.set('zmin', zr(1)).set('zmax', zr(2));
box.set('entitydim', 2).set('condition', condition);
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
