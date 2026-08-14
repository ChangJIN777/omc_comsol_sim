function [model,P] = buildBoomerangUnitCell(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.
%
% Builds the boomerang unit cell: the primitive rhombus of a hexagonal lattice
% (side a) with a three-pointed star hole - three arms, each w wide and r long,
% radiating from the cell centre at 120 deg spacing.
%
% OPTIONAL AIR DISK
%   P.add_airDisk     0/absent (default) = solid only, geometry unchanged.
%                     1 = add an air region on top of the slab, as in
%                     comsol_templates/BoomerangWithAirDisk.m.
%   P.airDiskHeight   Total height of that air region in metres, measured from
%                     the base of the remaining solid (z = 0 when P.mbevenz is
%                     nonzero, otherwise z = -th/2). P.airDiskH is accepted as
%                     an alias - see resolveAirDiskHeight below.
%   Only read when P.add_airDisk is on. See the CAUTION at the air disk block:
%   this file is on the mechanical path, and an air domain is an optical
%   construct.

%% read the input parameters
a = P.a;        % lattice constant; side of the rhombic cell
w = P.w;        % hole arm width (the narrowest etched feature)
th = P.th;      % full slab thickness in z
r = P.r;        % hole arm length, cell centre to arm tip
r1 = P.r1;      % fillet radius at the inner corners, where the arms meet
r2 = P.r2;      % fillet radius at the outer arm tips
abssym = abs(P.mbevenz);    % symmetry in the z direction

%% Create component 
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% create unit cell with boomerang geometry
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellWP.set('unite', true);

ucellplane = ucellWP.geom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
rec_1 = ucellWP.geom.feature.create('r1', 'Rectangle');
rec_1.set('pos', [a*(1/2+1/4) a*sqrt(3)/4+r/2]);
rec_1.set('base', 'center');
rec_1.set('size', [w r]);
rec_2 = ucellWP.geom.feature.create('r2', 'Rectangle');
rec_2.set('pos', [a*(1/2+1/4)-sqrt(3)*r/4 a*sqrt(3)/4-r/4]);
rec_2.set('base', 'center');
rec_2.set('rot', 120);
rec_2.set('size', [w r]);
rec_3 = ucellWP.geom.feature.create('r3', 'Rectangle');
rec_3.set('pos', [a*(1/2+1/4)+sqrt(3)*r/4 a*sqrt(3)/4-r/4]);
rec_3.set('base', 'center');
rec_3.set('rot', 240);
rec_3.set('size', [w r]);
composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'pol1-r1-r2-r3');

% add fillets 
hole_pos = [a*(1/2+1/4) a*sqrt(3)/4];
selection_width = 50e-9;
addFillet(P,ucellWP.geom,hole_pos,selection_width)
% ucellWP.geom.create('fil1', 'Fillet');
% ucellWP.geom.feature('fil1').set('radius', r1);
% ucellWP.geom.feature('fil1').selection('point').set('r1(1)', [3 4]);
% ucellWP.geom.feature('fil1').selection('point').set('r3(1)', [3 4]);
% ucellWP.geom.feature('fil1').selection('point').set('r2(1)', [3 4]);
% centerTriangle = ucellWP.geom.feature.create('pol2', 'Polygon');
% centerTriangle.set('source', 'table');
% centerTriangle.set('table', [a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); w+a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2+w/2 -(w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3)]);

extrude = ucellgeom.feature.create('ext1', 'Extrude');
extrude.setIndex('distance', th, 0);
extrude.selection('input').set({'wp1'});
ucellgeom.runAll;

%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('pol2', 'Polygon');
    symZPlane.set('source', 'table');
    symZPlane.set('table', [0 0; symW/2 (symW/2)*sqrt(3); symW*(3/2) (symW/2)*sqrt(3); symW 0; 0 0]);

    % extrude symmetry block
    ucellgeom.runCurrent;
    symZPlaneExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', symZth);

    % compose: unit cell - symmetry block
    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set('ext1');
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', ['ext1 - symZPlaneExt']);
    ucellgeom.runCurrent;


    % beam z-symmetry plane
    delta = 10e-9;
    ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
    % Unbounded in x and y. The rhombic cell spans x in [0, 3a/2] and
    % y in [0, a*sqrt(3)/2], so the previous +/-(a/2+delta) box sat mostly
    % outside the cell and clipped the z = 0 face it is meant to select.
    % Only the z extent has to be tight - that is what picks out the plane.
    ZsymSel.set('xmin', '-inf').set('xmax', 'inf');
    ZsymSel.set('ymin', '-inf').set('ymax', 'inf');
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end

%% Air disk (optional, P.add_airDisk)
% Adds an air region sitting on the slab, ported from
% comsol_templates/BoomerangWithAirDisk.m. That template builds it as a second
% work plane ('airDisk') carrying the SAME rhombic footprint as the slab base,
% extruded upward, and deliberately does NOT compose it with the solid: the
% default Form Union finalization partitions the overlap, so the etched tri-arm
% hole and the space above the slab both become air while the diamond keeps its
% own domain. Reusing the identical footprint also matters for the periodic BCs
% - an air block of any other outline would not tile with the Floquet pairs.
%
% Absent or zero P.add_airDisk leaves the geometry byte-for-byte as before, so
% every existing caller is unaffected.
%
% CAUTION, and please read before switching this on for a band-structure run.
% The only live caller of this file is runBands_2D (the MECHANICAL path); the
% optical boomerang runs go through buildBoomerangUnitCell_2D (2D and 3D) or
% buildBoomerangUnitCellStrip_v2 (1D) instead. runBands_2D hard-codes
% mbfem.b_domind = 1 and assigns both the material and the solid-mechanics
% selection to that one domain, so a second domain introduced here gets no
% material and no physics, and may renumber which domain is the diamond. The
% added boundaries can also change what the bndindex lookups below return for
% the Floquet and symmetry selections. An air disk is an optical construct; if
% that is what you are after, buildBoomerangUnitCell_2D is probably the file you
% want. It is implemented here because it was asked for, defaulted off so it
% cannot surprise a mechanical study.
if isfield(P,'add_airDisk') && P.add_airDisk
    airDiskHeight = resolveAirDiskHeight(P);

    % Base of the air region. The slab is extruded from -th/2 through +th/2,
    % but the z-symmetry block above removes everything below z = 0 when
    % P.mbevenz is nonzero - so the air has to start from whichever face is
    % actually the bottom of the remaining solid, or it would begin inside it.
    if abs(P.mbevenz)
        zAirBase = 0;
    else
        zAirBase = -th/2;
    end

    % airDiskHeight is the TOTAL height of the air region measured from that
    % base, matching the template's flat 5000[um] extrude from z = 0 - not the
    % clearance above the top face. Warn rather than error if it fails to clear
    % the solid: the geometry still builds, it just contains no air above the
    % slab, which is a silent way to get meaningless optical results.
    if zAirBase + airDiskHeight <= th/2
        warning('buildBoomerangUnitCell:airDiskTooShort', ...
            ['Air disk height %g nm measured from z = %g nm does not reach ' ...
             'above the slab top at %g nm, so no air lies above the solid.'], ...
            airDiskHeight*1e9, zAirBase*1e9, (th/2)*1e9);
    end

    airWP = ucellgeom.feature.create('airDiskWP', 'WorkPlane');
    airWP.label('airDisk');
    airWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', zAirBase);
    airPlane = airWP.geom.feature.create('pol_air', 'Polygon');
    airPlane.set('source', 'table');
    airPlane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);

    airExt = ucellgeom.feature.create('airDiskExt', 'Extrude');
    airExt.setIndex('distance', airDiskHeight, 0);
    % Input selection set explicitly. The symmetry-block extrude above omits
    % it and relies on the default, which silently picks up the wrong work
    % plane once there is more than one candidate in the sequence.
    airExt.selection('input').set({'airDiskWP'});
    ucellgeom.runCurrent;
end

%% Periodic-boundary selections (named)
% Cumulative selections, so the Floquet BCs can be picked by name -
% 'geom1_xboundaries_bnd' and 'geom1_yboundaries_bnd' - the way runBands.m
% does it, instead of by the P.xEnd*/P.yEnd* indices below. Index lookups
% renumber whenever a geometry feature is added; see the air disk CAUTION.
%
% This cell is the primitive rhombus of a hexagonal lattice, not a rectangle,
% so the pair separated by a1 = (a,0) - the "x" pair - are the two 60 deg
% SLANTED faces. A thin axis-aligned box cannot isolate a slanted face, so
% each is caught by a box spanning its full x extent (left face x in
% [0, a/2], right face x in [a, 3a/2]) with condition 'inside'. The a2 pair
% (bottom y = 0, top y = a*sqrt(3)/2) is axis-aligned, so a thin y slab does.
sel_delta = 10e-9;

x_boundary_boxsel_l = ucellgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_boxsel_l.set('entitydim', 2);
x_boundary_boxsel_l.set('xmin', -sel_delta).set('xmax', a/2+sel_delta);
x_boundary_boxsel_l.set('inputent', 'all');
x_boundary_boxsel_l.set('condition', 'inside');

x_boundary_boxsel_r = ucellgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_boxsel_r.set('entitydim', 2);
x_boundary_boxsel_r.set('xmin', a-sel_delta).set('xmax', 3*a/2+sel_delta);
x_boundary_boxsel_r.set('inputent', 'all');
x_boundary_boxsel_r.set('condition', 'inside');

ucellgeom.selection.create('xboundaries','CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
x_boundary_boxsel_l.set('contributeto','xboundaries');
x_boundary_boxsel_r.set('contributeto','xboundaries');

y_boundary_boxsel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
y_boundary_boxsel_b.set('entitydim', 2);
y_boundary_boxsel_b.set('ymin', -sel_delta).set('ymax', sel_delta);
y_boundary_boxsel_b.set('inputent', 'all');
y_boundary_boxsel_b.set('condition', 'inside');

y_boundary_boxsel_t = ucellgeom.feature.create('y_boundary_boxsel_t', 'BoxSelection');
y_boundary_boxsel_t.set('entitydim', 2);
y_boundary_boxsel_t.set('ymin', (a/2)*sqrt(3)-sel_delta).set('ymax', (a/2)*sqrt(3)+sel_delta);
y_boundary_boxsel_t.set('inputent', 'all');
y_boundary_boxsel_t.set('condition', 'inside');

ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
y_boundary_boxsel_b.set('contributeto','yboundaries');
y_boundary_boxsel_t.set('contributeto','yboundaries');

ucellgeom.runAll;

%% Making selections (manual)
mphgeom(model);
P.xEnd1 =  bndindex(ucellgeom, [0 0 0], [sqrt(3)*a/2 -a/2 0]);
P.xEnd2 = bndindex(ucellgeom, [a 0 0], [sqrt(3)*a/2 -a/2 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 0 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [a (a/2)*sqrt(3) 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);

% Cross-check the named selections against the index lookups above. The x
% boxes span half the cell each, so if a hole arm reaches into one of those
% bands (large r/a) the 'inside' condition would swallow an interior face
% too, and a Floquet condition would be applied to the wrong surface. Warn
% rather than error - the geometry is still valid, the selection is not.
checkNamedBndSel(model, ucellname, 'xboundaries', [P.xEnd1 P.xEnd2]);
checkNamedBndSel(model, ucellname, 'yboundaries', [P.yEnd1 P.yEnd2]);
% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end

function airDiskHeight = resolveAirDiskHeight(P)
%RESOLVEAIRDISKHEIGHT Air disk height [m], from P.airDiskHeight or P.airDiskH.
%
% P.airDiskHeight is the documented name for this builder. P.airDiskH is
% accepted as an alias because it is what the rest of this directory already
% uses - buildHoleUnitCell.m, buildHoleStrip_3D.m and
% buildBoomerangUnitCellStrip_v2.m all read P.airDiskH, and test_BoomerangStrip.m,
% test_Hole_Strip.m, test_hole_unitCell.m, optimize_hole_unitCell.m and
% boomerang_optimize_sweep_diamond.m all set it. Rejecting it would mean those
% configs set an air disk height that this file silently ignored. When both are
% present P.airDiskHeight wins, since it is the more specific request.

if isfield(P,'airDiskHeight')
    airDiskHeight = P.airDiskHeight;
elseif isfield(P,'airDiskH')
    airDiskHeight = P.airDiskH;
else
    error('buildBoomerangUnitCell:noAirDiskHeight', ...
        ['P.add_airDisk is enabled but no height was given. Set ' ...
         'P.airDiskHeight (or P.airDiskH) to the air region height in metres.']);
end

% Caught here rather than inside COMSOL, where a zero or negative extrude
% distance surfaces as an opaque Java geometry error several calls deeper.
validateattributes(airDiskHeight, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'air disk height');
end

% -------------------------------------------------------------------------

function checkNamedBndSel(model, ucellname, selname, expected)
%CHECKNAMEDBNDSEL Warn if a cumulative selection misses the expected faces.
%
% selname   geometry cumulative selection, e.g. 'xboundaries'
% expected  boundary indices the bndindex lookups found for that pair
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
        warning('buildBoomerangUnitCell:selectionUnreadable', ...
            ['Could not read selection %s (%s), so it was not verified ' ...
             'against the bndindex lookup.'], tag, err.message);
        return
    end
end

got = sort(got(:))';
expected = sort(double(expected(:)))';

if isempty(got)
    warning('buildBoomerangUnitCell:selectionEmpty', ...
        ['%s resolved to no boundaries. The Floquet condition using it ' ...
         'would be applied to nothing.'], tag);
elseif ~isequal(got, expected)
    warning('buildBoomerangUnitCell:selectionMismatch', ...
        ['%s resolved to boundaries [%s] but bndindex gives [%s]. Check ' ...
         'the geometry before trusting the periodic BCs - a hole arm may ' ...
         'reach into the selection box.'], ...
        tag, num2str(got), num2str(expected));
end
end

% -------------------------------------------------------------------------

function addFillet(P,ucellgeom,hole_pos,selection_width)
    w = P.w;
    r = P.r;
    r1 = P.r1;
    r2 = P.r2;
    disksel1_label = sprintf('h_disksel1');
    disksel2_label = sprintf('h_disksel2');
    fillet1_label = sprintf('h_fil1');
    fillet2_label = sprintf('h_fil2');
    % hole 
    disksel1 = ucellgeom.feature.create(disksel1_label, 'DiskSelection');
    disksel1.set('entitydim', 0);
    disksel1.set('posx', hole_pos(1));
    disksel1.set('posy', hole_pos(2));
    disksel1.set('r', w);
    disksel1.set('rin', w/(2*sqrt(2)));
    disksel1.set('condition', 'allvertices');
    fil1 = ucellgeom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
    disksel2 = ucellgeom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', hole_pos(2));
    disksel2.set('posx', hole_pos(1));
    disksel2.set('r', r+selection_width/2);
    disksel2.set('rin', r-selection_width/2);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellgeom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end


