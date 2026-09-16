function [model,P] = BuildNanobeamCrossShieldFEM(model,P)
%BUILDNANOBEAMCROSSSHIELDFEM  Nanobeam cavity + 2D cross-cell phononic shield + wrap-around mechanical PML.
%
%   [model,P] = BuildNanobeamCrossShieldFEM(model,P)
%
%   Geometry builder for RunNanobeamFEM, selected by P.celltype = 'crossShield'.
%   It is BuildNanobeamFEM.m with the +x beam end terminated by a two
%   dimensional array of cross unit cells instead of by a fixed constraint or
%   by a single PML end cap, so that the radiative mechanical Q of the cavity
%   mode can be extracted from the complex eigenfrequency.
%
%   WHAT IT BUILDS  (symmetric cavity: quarter model in x,y; eighth with z)
%
%       x = 0 ................ cavity centre, mirror plane (P.mevenx)
%       x in [0, xB1] ........ half beam with elliptical holes, y in [0, w/2]
%       x in [xB1, xS0] ...... optional homogeneous clamp pad, y in [0, yS1]
%       x in [xS0, xS1] ...... cross shield, nShieldX by nShieldY cells
%       x in [xS1, xP1] ...... mechanical PML, +x arm
%       y in [yS1, yP1] ...... mechanical PML, +y arm    (x in [xB1, xS1])
%       corner ............... mechanical PML, +x/+y     (2-direction stretch)
%
%   with
%       xB1 = P.beamLenHalf                 (P.beamLenPhMHalf if nphonmir > 0)
%       xS0 = xB1 + P.shieldPadLen
%       xS1 = xS0 + P.nShieldX*P.aShield
%       yS1 =       P.nShieldY*P.aShield
%       xP1 = xS1 + P.PMLLen,   yP1 = yS1 + P.PMLLenY
%
%   CROSS CELL GEOMETRY
%     Ported from buildCrossStrip.m, which must NOT be called directly: it
%     reads P.a, P.h, P.w, P.r1, P.r2, and in this pipeline P.a is the HOLE
%     lattice constant and P.w is the BEAM WIDTH. It also creates its own
%     'comp1'/'geom1'. The loop below is the same construction with explicit
%     numeric arguments and a 2D (i,j) index.
%
%     Cross void centres, square lattice, pitch P.aShield in BOTH x and y:
%
%         x_i = xS0 + (i - 1/2)*aShield,   i = 1 .. P.nShieldX
%         y_j =       (j - 1/2)*aShield,   j = 1 .. P.nShieldY
%
%     Each void is the union of two crossed rectangles, [hShield wShield] and
%     [wShield hShield], subtracted from the slab footprint.
%
%     The resulting solid is a square lattice of NODES of side
%     (aShield - wShield), centred on the CELL CORNERS, joined by LIGAMENTS of
%     width (aShield - hShield) and length wShield.
%
%   WHY THIS CELL PHASING
%     The void never reaches a cell edge (its extreme reach is hShield/2 <
%     aShield/2 in both directions), so every cell-edge plane is entirely
%     solid. Three things follow, and all three are load bearing:
%
%       1. y = 0 is a full-solid mirror plane bisecting a row of nodes, so the
%          half-beam y symmetry survives the shield unchanged.
%       2. x = xS0 is fully solid, so the beam end face unions cleanly onto
%          the shield with no dangling or non-manifold geometry.
%       3. x = xS1 and y = yS1, the shield/PML interfaces, are clean planar
%          rectangles with no void opening onto them. COMSOL's Cartesian PML
%          infers its stretching direction from the PML domain's INNER
%          boundary; a ragged, void-punctured interface is where that
%          inference goes wrong.
%
%   P.nShieldY IS A HALF-MODEL COUNT
%     P.nShieldY is the number of FULL cells in the simulated y >= 0 half. The
%     physical device has 2*P.nShieldY rows, and the mirror plane bisects a row
%     of nodes. P.shield.nYPhys = 2*P.nShieldY is written back into P and
%     appears in the filename (nsyH<half>_nsyP<physical>) so it cannot be
%     misread later.
%
%     The alternative convention, a half cell straddling y = 0, is NOT
%     supported and is not a bug: it would put the beam axis on a VOID centre
%     line, so the y = 0 plane through the first cell would be solid only over
%     |y| < (aShield-hShield)/2 and the beam would butt into a membrane of that
%     thickness.
%
%   THE CLAMP COMPROMISE - READ THIS BEFORE TRUSTING A Q
%     With P.shieldPadLen = 0 (the default, lattice faithful) the beam end face
%     lands on a half node of half-width (aShield - wShield)/2, and on the
%     half-wall fin of thickness (aShield - hShield)/2 for everything beyond
%     that. For the legacy cross parameters (a=894, h=863, w=416 nm) against a
%     750 nm beam that is 36% of the beam end half-width landing on a 15.5 nm
%     wall. auditCrossShield below computes the actual fraction and WARNS. It
%     never changes the geometry.
%
%     The fix, when you want it, is P.shieldPadLen > 0: a homogeneous clamp
%     slab spanning the full shield height, which is what a real device has.
%     Recommended production value is about aShield/2. It is implemented and
%     tested here; it is simply not the default.
%
%   FILLETS
%     P.shieldFillet defaults to 'none' and emits no DiskSelection/Fillet
%     features at all. Two other modes exist:
%
%       'legacy'   - byte for byte the annuli of buildCrossUnitCell.m:129-148,
%                    disksel1 = [w/2, 1.5*w] and disksel2 = [h/2-5nm, h/2+5nm].
%                    EXPECT THIS TO FAIL FOR nx,ny > 1. In a 2D array each
%                    cell's disksel1 captures 8 FOREIGN vertices (two from each
%                    of the four orthogonal neighbours, at radius
%                    sqrt((a-h/2)^2 + (w/2)^2) = 507 nm for the legacy
%                    parameters, inside the [208, 624] nm annulus) - exactly
%                    twice the bleed buildCrossStrip.m:119-136 documents for a
%                    1D strip. Fillets run in sequence, so a cell rounds its
%                    neighbours' tips before those neighbours' own selections
%                    are evaluated, and the neighbour then tries to fillet an
%                    already-filleted corner.
%       'analytic' - annuli centred on the vertices the names claim: the four
%                    inner concave corners at w/sqrt(2) get r1, the eight
%                    arm-tip corners at sqrt(h^2+w^2)/2 get r2. Clears the
%                    nearest foreign vertex by 23 nm for the legacy parameters;
%                    the margin is computed and warned on.
%
%     'analytic' CHANGES THE GEOMETRY even at identical r1/r2, because legacy
%     rounds all 12 corners with r1 and none with r2. It is not a drop-in match
%     for already-simulated or fabricated cross designs.
%
%   REQUIRED P FIELDS  (SI, metres - see CLAUDE.md)
%     Everything BuildNanobeamFEM needs, plus:
%       P.celltype     'crossShield'
%       P.aShield      cross lattice pitch, both x and y [m]
%       P.hShield      cross arm length [m]  (hard error if >= aShield)
%       P.wShield      cross arm width  [m]
%       P.nShieldX     cells along x (integer >= 1)
%       P.nShieldY     cells along y in the SIMULATED HALF (integer >= 1)
%
%   OPTIONAL P FIELDS
%       P.shieldPadLen    clamp slab length [m].                 Default 0
%       P.shieldFillet    'none' | 'analytic' | 'legacy'.        Default 'none'
%       P.r1Shield        inner-corner fillet radius [m].        Default 0
%       P.r2Shield        arm-tip fillet radius [m].             Default 0
%       P.shieldMinFeature  lithography limit for the audit [m]. Default 50e-9
%       P.shieldHmax      max mesh element in the shield [m].
%                         Default min((a-h)/3, th/3) - three elements across a
%                         ligament. Consumed by SolveNanobeamFEM.
%       P.PMLLenY         +y PML thickness [m].               Default P.PMLLen
%
%   PROVIDED TO THE CALLER  (the contract SetupNanobeamFEM/SolveNanobeamFEM use)
%       P.domSel.beam     beam + pad + shield  (material and smech selection)
%       P.domSel.shield   shield only          (mesh sizing)
%       P.domSel.PML      the three PML frame domains
%       P.bndSel.beamXsym x = 0  cavity mirror
%       P.bndSel.beamYsym y = 0  faces of beam + pad + shield
%       P.bndSel.PMLYsym  y = 0  face of the +x PML arm
%       P.bndSel.beamZsym z = 0  faces of beam + pad + shield   (if symZOn)
%       P.bndSel.PMLZsym  z = 0  faces of the PML frame         (if symZOn)
%       P.bndSel.beamXend outer termination faces - the shield perimeter when
%                         the PML is off, the PML outer faces when it is on.
%                         SetupNanobeamFEM.m:221-223 only reads it in the
%                         no-PML branch, where clamping the shield perimeter is
%                         the right thing.
%       P.bndSel.shieldXend, P.bndSel.shieldYend, P.bndSel.PMLcurv  (debug)
%       P.shield          struct of computed extents and cell centres
%
%   NOT SUPPORTED IN v1 (all hard errors, see the checks below)
%       P.asymCav      - needs a mirrored shield and three more PML blocks at
%                        -x. addCrossShield takes an (xEdge, sgn) pair so the
%                        left end is a sign flip rather than a rewrite.
%       P.xsect ~= 'rect' - 'tri'/'isoFit' force P.mevenz = 0 and make the beam
%                        an intersection with a prism; the junction to a
%                        prismatic cross slab is undefined.
%       P.solveOpt     - the air cylinder (radius ~ 2*lambda + w/2) would
%                        swallow the shield and the PML frame.
%       P.meveny == 0  - there is no full-y build path in this pipeline.
%
%   KNOWN BUG IN THE LEGACY PML PATH, DELIBERATELY NOT TOUCHED HERE
%     BuildNanobeamFEM.m:508, 518 and 538 build PMLtopSel / PMLbotSel /
%     PMLxoutSel with ymax = max(w)/2 + delta and condition 'allvertices', but
%     the PML pad created at line 241 extends to y = P.PMLLen, which is orders
%     of magnitude larger. Those selections therefore resolve to EMPTY, so with
%     P.mevenz = 1 the z symmetry BC is applied to nothing on the PML and its
%     z = 0 face is free instead of mirrored. That file is left exactly as it
%     is - fixing it would change results for every existing caller. This
%     builder computes yMaxAll from the real footprint instead, and
%     assertNonEmptySel warns on every selection it emits so the same class of
%     failure cannot pass silently again.
%
%   See also BUILDNANOBEAMFEM, BUILDCROSSSTRIP, BUILDCROSSUNITCELL,
%            ISCROSSFABRICABLE, RUNNANOBEAMFEM, SETUPNANOBEAMFEM.

%% Read, default and validate the parameters before COMSOL sees anything
% A bad size or a zero extrude distance surfaces from COMSOL as an opaque Java
% geometry error several calls deeper (buildCrossStrip.m:362 makes the same
% argument).
P = readCrossShieldParams(P);

wid = P.w;
thi = P.th;

geom    = P.geomHalf;
beamLen = P.beamLenHalf;
P.xc    = 0;

% With a phononic mirror appended, CreateNanobeamGeom.m:424 derives
% beamLenHalf from index nholes+nwgm, which is not the last hole, so the outer
% mirror holes would fall outside the beam rectangle and be silently dropped.
if isfield(P,'nphonmir') && P.nphonmir > 0
    if ~isfield(P,'beamLenPhMHalf')
        error('BuildNanobeamCrossShieldFEM:noPhMLength', ...
            'P.nphonmir = %d but P.beamLenPhMHalf is missing.', P.nphonmir);
    end
    if P.beamLenPhMHalf > beamLen + 1e-15
        error('BuildNanobeamCrossShieldFEM:beamShorterThanHoles', ...
            ['P.nphonmir = %d puts holes out to x = %.1f nm but ' ...
             'P.beamLenHalf is only %.1f nm, so the outer phononic-mirror ' ...
             'holes fall outside the beam. Set P.nphonmir = 0, or use a ' ...
             'geometry generator that reports the full length.'], ...
            P.nphonmir, P.beamLenPhMHalf*1e9, beamLen*1e9);
    end
end

nholes = size(geom,1);
hx     = geom(:,1)';
hy     = geom(:,2)';
xpos   = geom(:,3)';
ypos   = geom(:,4)';

%% Derived extents - these replace BuildNanobeamFEM's totLen/maxWid/xL
% BuildNanobeamFEM tracks the model bounding box with three running scalars
% (xL, totLen, maxWid) that are accumulated as features are added and then
% reused as SIZES in the z-symmetry block. That bookkeeping is what breaks once
% the structure is wider than the beam: see the KNOWN BUG note above. Compute
% the box explicitly instead, once, from the geometry that is about to be
% built.
as = P.aShield;
hs = P.hShield;
ws = P.wShield;
nx = P.nShieldX;
ny = P.nShieldY;

xB1 = beamLen;                      % beam end / clamp pad start
xS0 = xB1 + P.shieldPadLen;         % shield left edge
xS1 = xS0 + nx*as;                  % shield right edge
yS1 = ny*as;                        % shield top edge

usePML = P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML;
if usePML
    xP1 = xS1 + P.PMLLen;
    yP1 = yS1 + P.PMLLenY;
else
    xP1 = xS1;
    yP1 = yS1;
end

xMinAll = 0;
xMaxAll = xP1;
yMaxAll = yP1;

xcentres = xS0 + ((1:nx) - 0.5)*as;
ycentres =       ((1:ny) - 0.5)*as;

% z symmetry is only meaningful for the rectangular cross-section, and v1 is
% rect only, so the guard collapses to the mechanical condition.
symZOn = abs(P.mevenz) > 0;
if symZOn
    zLoEff = 0;
else
    zLoEff = -thi/2;
end

P.shield = struct( ...
    'a',            as, ...
    'h',            hs, ...
    'w',            ws, ...
    'nX',           nx, ...
    'nY',           ny, ...
    'nYPhys',       2*ny, ...
    'padLen',       P.shieldPadLen, ...
    'xB1',          xB1, 'xS0', xS0, 'xS1', xS1, 'yS1', yS1, ...
    'xP1',          xP1, 'yP1', yP1, ...
    'xcentres',     xcentres, ...
    'ycentres',     ycentres, ...
    'nodeWidth',    as - ws, ...
    'wallWidth',    as - hs, ...
    'solidFrac',    1 - (2*hs*ws - ws^2)/as^2);

%% Pure-MATLAB audit - warn only, never alters the geometry
auditCrossShield(P, wid, thi, usePML);

%% Create component
comp = model.modelNode.create('comp');
comp.label('Nanobeam FEM simulation');

P.geomname = 'beam';
beamgeom = model.geom.create(P.geomname, 3);
beamgeom.label('Nanobeam with cross phononic shield');

%% Beam with rectangular cross-section and elliptical holes
beamWP = beamgeom.feature.create('beamWP', 'WorkPlane');
beamWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -thi/2);

beamplane = beamWP.geom.feature.create('beamplane', 'Rectangle');
beamplane.set('type', 'solid').set('base', 'corner');
beamplane.set('pos', [0 0]).set('size', [beamLen wid/2]);
beamgeom.runCurrent;

holeList = cell(1, nholes);
for k = 1:nholes
    holeID = ['hole_', num2str(k)];
    hole = beamWP.geom.feature.create(holeID, 'Ellipse');
    hole.set('pos', [xpos(k) ypos(k)]);
    hole.set('semiaxes', [hx(k)/2 hy(k)/2]);
    holeList{k} = holeID;
end
beamgeom.runCurrent;

beamComp = beamWP.geom.feature.create('beamComp', 'Compose');
beamComp.selection('input').set([{'beamplane'}, holeList]);
beamComp.set('formula', ['beamplane - ', strjoin(holeList, ' - ')]);
beamComp.set('intbnd', false);
beamgeom.runCurrent;

% Extrude with an EXPLICIT input selection. BuildNanobeamFEM.m:125-127 omits
% this and relies on the default picking the only work plane; buildCrossStrip.m
% :233-236 documents why that stops working - "the default input silently picks
% the wrong work plane once the sequence holds more than one candidate". This
% file adds a second work plane, so both extrudes must set it.
beamHoles = beamgeom.feature.create('beamHoles', 'Extrude');
beamHoles.set('distance', thi);
beamHoles.selection('input').set({'beamWP'});
beamgeom.runCurrent;

finBeamTag = 'beamHoles';

%% Cross shield at the +x end
% sgn = +1 builds to the right of xEdge. A mirrored left-end shield for the
% asymmetric cavity is the sgn = -1 call; it is not wired up in v1 (P.asymCav
% is rejected in readCrossShieldParams) but the emitter takes the argument so
% adding it later is a sign flip, not a rewrite.
shieldTag = addCrossShield(beamgeom, P, xB1, +1, thi);

% intbnd = true (KEEP interior boundaries) is set explicitly, not left to the
% default, because the whole downstream contract depends on it: the beam and
% the shield must remain SEPARATE domains so that P.domSel.shield can carry its
% own mesh Size node. Continuity across an interior boundary is automatic
% within one physics, so nothing is lost mechanically.
beamShield = beamgeom.feature.create('beamShield', 'Compose');
beamShield.selection('input').set(finBeamTag);
beamShield.selection('input').set(shieldTag);
beamShield.set('formula', [finBeamTag, ' + ', shieldTag]);
beamShield.set('intbnd', true);
beamgeom.runCurrent;
finBeamTag = 'beamShield';

%% Mechanical PML - L frame wrapping the shield perimeter
% Which sides: the structure is a suspended slab in vacuum, so mechanical
% energy leaks only IN PLANE, through material. There is no material above or
% below, hence no +/-z PML; x = 0 and y = 0 are mirror planes, hence no -x/-y
% PML; and the beam's own +y flank for x < xS0 is free vacuum, not a radiating
% interface. That leaves +x, +y and the corner.
%
% Three blocks, NOT two, and interior boundaries kept, so that the PML resolves
% to THREE domains. COMSOL's box PML assigns stretching per domain from that
% domain's inner boundary: x only for the +x arm, y only for the +y arm, and x
% AND y for the corner. Fold the corner into an arm and it stretches in one
% direction only, and leaks.
displayPMLStr = '';
if usePML
    zb = -thi/2;

    PMLxArm = beamgeom.feature.create('PMLxArm', 'Block');
    PMLxArm.set('base', 'corner');
    PMLxArm.set('pos',  [xS1, 0, zb]);
    PMLxArm.set('size', [P.PMLLen, yS1, thi]);

    % The +y arm spans the whole slab footprint in x, including the clamp pad,
    % because the pad's +y face radiates too.
    PMLyArm = beamgeom.feature.create('PMLyArm', 'Block');
    PMLyArm.set('base', 'corner');
    PMLyArm.set('pos',  [xB1, yS1, zb]);
    PMLyArm.set('size', [xS1 - xB1, P.PMLLenY, thi]);

    PMLcorner = beamgeom.feature.create('PMLcorner', 'Block');
    PMLcorner.set('base', 'corner');
    PMLcorner.set('pos',  [xS1, yS1, zb]);
    PMLcorner.set('size', [P.PMLLen, P.PMLLenY, thi]);

    beamgeom.runCurrent;

    % intbnd = true here is what keeps the corner a domain of its own, which is
    % the whole reason for building three blocks instead of two. The
    % assertNonEmptySel audit below checks that the count really is 3.
    PMLAll = beamgeom.feature.create('MechPML', 'Compose');
    PMLAll.selection('input').set({'PMLxArm','PMLyArm','PMLcorner'});
    PMLAll.set('formula', 'PMLxArm + PMLyArm + PMLcorner');
    PMLAll.set('intbnd', true);
    beamgeom.runCurrent;

    beamPML = beamgeom.feature.create('beamShieldPML', 'Compose');
    beamPML.selection('input').set(finBeamTag);
    beamPML.selection('input').set('MechPML');
    beamPML.set('formula', [finBeamTag, ' + MechPML']);
    beamPML.set('intbnd', true);
    beamgeom.runCurrent;

    finBeamTag = 'beamShieldPML';
    displayPMLStr = ', mech PML (L frame around shield)';
end

%% Symmetry in z
% Footprint from the computed box, not from accumulated totLen/maxWid scalars.
if symZOn
    symZWP = beamgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -thi/2);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [xMinAll 0]).set('size', [xMaxAll - xMinAll, yMaxAll]);
    beamgeom.runCurrent;

    symZPlaneExt = beamgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', thi/2);
    symZPlaneExt.selection('input').set({'symZWP'});

    symZComp = beamgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set(finBeamTag);
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', [finBeamTag, ' - symZPlaneExt']);
    beamgeom.runCurrent;
end

beamgeom.run;

disp(['Geometry created - nanobeam, rectangular cross-section, ', ...
      num2str(nx), 'x', num2str(ny), ' cross shield (', ...
      num2str(2*ny), ' physical rows)', displayPMLStr]);

%% Domain selections
d = 10e-9;

beamSel = beamgeom.create('beamSel', 'BoxSelection');
beamSel.set('xmin', -d).set('xmax', xB1 + d);
beamSel.set('ymin', -d).set('ymax', wid/2 + d);
beamSel.set('zmin', zLoEff - d).set('zmax', thi/2 + d);
beamSel.set('entitydim', 3).set('condition', 'inside');
beamgeom.runCurrent;

shieldSel = beamgeom.create('shieldSel', 'BoxSelection');
shieldSel.set('xmin', xS0 - d).set('xmax', xS1 + d);
shieldSel.set('ymin', -d).set('ymax', yS1 + d);
shieldSel.set('zmin', zLoEff - d).set('zmax', thi/2 + d);
shieldSel.set('entitydim', 3).set('condition', 'inside');
beamgeom.runCurrent;
P.domSel.shield = readSel(model, P.geomname, 'shieldSel');

% beam + pad + shield. This is what gets the material and the Solid Mechanics
% physics; BuildNanobeamFEM's beamSel stops at ymax = max(w)/2, which would
% leave the shield with no material at all.
structSel = beamgeom.create('structSel', 'BoxSelection');
structSel.set('xmin', -d).set('xmax', xS1 + d);
structSel.set('ymin', -d).set('ymax', yS1 + d);
structSel.set('zmin', zLoEff - d).set('zmax', thi/2 + d);
structSel.set('entitydim', 3).set('condition', 'inside');
beamgeom.runCurrent;
P.domSel.beam = readSel(model, P.geomname, 'structSel');

if usePML
    % PML domains are selected with 'intersects' PROBE boxes, one per arm.
    %
    % Why probes rather than one frame-shaped 'intersects' box: the PML is an
    % L (+x arm, +y arm, corner), which no single box bounds without also
    % containing the shield. An 'intersects' box drawn around the whole frame
    % would ALSO select the shield and the beam, because a domain that merely
    % shares the x = xS1 or y = yS1 face with the box still intersects it.
    % A small box placed strictly inside one arm intersects that arm and
    % nothing else, and unlike the previous 'inside' + DifferenceSelection
    % construction it does not depend on the outer PML bounds (xP1, yP1,
    % zLoEff) being exactly right.
    %
    % Arm extents come from the Block features above:
    %   PMLxArm    x in [xS1, xP1],  y in [0,   yS1]
    %   PMLyArm    x in [xB1, xS1],  y in [yS1, yP1]
    %   PMLcorner  x in [xS1, xP1],  y in [yS1, yP1]
    % Each probe sits at the centre of its arm, so the margin to the nearest
    % arm face is half the smallest arm dimension - microns against d = 10 nm.
    zMid = 0.5*(zLoEff + thi/2);
    probes = { ...
        'PMLxArmSel', 0.5*(xS1 + xP1), 0.5*yS1,         zMid; ...
        'PMLyArmSel', 0.5*(xB1 + xS1), 0.5*(yS1 + yP1), zMid; ...
        'PMLcornSel', 0.5*(xS1 + xP1), 0.5*(yS1 + yP1), zMid};

    PMLinds = [];
    for k = 1:size(probes,1)
        tag = probes{k,1};
        pk  = beamgeom.create(tag, 'BoxSelection');
        pk.set('xmin', probes{k,2} - d).set('xmax', probes{k,2} + d);
        pk.set('ymin', probes{k,3} - d).set('ymax', probes{k,3} + d);
        pk.set('zmin', probes{k,4} - d).set('zmax', probes{k,4} + d);
        pk.set('entitydim', 3).set('condition', 'intersects');
        beamgeom.runCurrent;
        ik = readSel(model, P.geomname, tag);
        assertNonEmptySel(['PML probe ', tag], ik);
        PMLinds = [PMLinds, ik];                                    %#ok<AGROW>
    end
    P.domSel.PML = unique(PMLinds);
    disp(['PML domain indices: ', num2str(P.domSel.PML)]);

    % Warn-only cross-check against the 'inside' + DifferenceSelection
    % construction this replaced, so a disagreement surfaces on the first run
    % instead of as an implausible Q. Costs one extra selection feature.
    allSel = beamgeom.create('allSel', 'BoxSelection');
    allSel.set('xmin', -d).set('xmax', xP1 + d);
    allSel.set('ymin', -d).set('ymax', yP1 + d);
    allSel.set('zmin', zLoEff - d).set('zmax', thi/2 + d);
    allSel.set('entitydim', 3).set('condition', 'inside');
    beamgeom.runCurrent;

    PMLSel = beamgeom.create('PMLSel', 'DifferenceSelection');
    PMLSel.set('entitydim', 3).set('add', 'allSel').set('subtract', 'structSel');
    beamgeom.runCurrent;
    PMLindsInside = readSel(model, P.geomname, 'PMLSel');
    if ~isequal(sort(P.domSel.PML), sort(PMLindsInside))
        warning('BuildNanobeamCrossShieldFEM:PMLSelMismatch', ...
            ['The intersects probes selected PML domains [%s] but the ' ...
             'inside/difference construction gives [%s]. Inspect the ' ...
             'geometry in the GUI before trusting the PML.'], ...
            num2str(P.domSel.PML), num2str(PMLindsInside));
    end
end

%% Boundary selections
% Flat planes: box selection with condition that ALL VERTICES are in the box.
% BuildNanobeamFEM.m:422-429 uses a 5 nm 'intersects' probe for the beam z
% symmetry plane, which by construction can only ever find the beam's own
% faces - not the shield's, and not the pad's. Use full-footprint boxes here.
d2 = 5e-9;

beamXsymSel = beamgeom.create('beamXsymSel', 'BoxSelection');
beamXsymSel.set('xmin', -d2).set('xmax', d2);
beamXsymSel.set('ymin', -d2).set('ymax', wid/2 + d2);
beamXsymSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
beamXsymSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
P.bndSel.beamXsym = readSel(model, P.geomname, 'beamXsymSel');

% y = 0 faces of beam + pad + shield. Deliberately stops at xS1 so it does not
% also claim the PML's y = 0 face: SetupNanobeamFEM.m:246-249 concatenates this
% with P.bndSel.PMLYsym into a single selection.set call.
structYsymSel = beamgeom.create('structYsymSel', 'BoxSelection');
structYsymSel.set('xmin', -d2).set('xmax', xS1 + d2);
structYsymSel.set('ymin', -d2).set('ymax', d2);
structYsymSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
structYsymSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
P.bndSel.beamYsym = readSel(model, P.geomname, 'structYsymSel');

if symZOn
    structZsymSel = beamgeom.create('structZsymSel', 'BoxSelection');
    structZsymSel.set('xmin', -d2).set('xmax', xS1 + d2);
    structZsymSel.set('ymin', -d2).set('ymax', yS1 + d2);
    structZsymSel.set('zmin', -d2).set('zmax', d2);
    structZsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    P.bndSel.beamZsym = readSel(model, P.geomname, 'structZsymSel');
end

% Shield perimeter - the outer termination when there is no PML, and a useful
% debug handle when there is.
shieldXendSel = beamgeom.create('shieldXendSel', 'BoxSelection');
shieldXendSel.set('xmin', xS1 - d2).set('xmax', xS1 + d2);
shieldXendSel.set('ymin', -d2).set('ymax', yS1 + d2);
shieldXendSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
shieldXendSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
P.bndSel.shieldXend = readSel(model, P.geomname, 'shieldXendSel');

shieldYendSel = beamgeom.create('shieldYendSel', 'BoxSelection');
shieldYendSel.set('xmin', xB1 - d2).set('xmax', xS1 + d2);
shieldYendSel.set('ymin', yS1 - d2).set('ymax', yS1 + d2);
shieldYendSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
shieldYendSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
P.bndSel.shieldYend = readSel(model, P.geomname, 'shieldYendSel');

if usePML
    PMLYsymSel = beamgeom.create('PMLYsymSel', 'BoxSelection');
    PMLYsymSel.set('xmin', xS1 - d2).set('xmax', xP1 + d2);
    PMLYsymSel.set('ymin', -d2).set('ymax', d2);
    PMLYsymSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
    PMLYsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    P.bndSel.PMLYsym = readSel(model, P.geomname, 'PMLYsymSel');

    if symZOn
        allZsymSel = beamgeom.create('allZsymSel', 'BoxSelection');
        allZsymSel.set('xmin', -d2).set('xmax', xP1 + d2);
        allZsymSel.set('ymin', -d2).set('ymax', yP1 + d2);
        allZsymSel.set('zmin', -d2).set('zmax', d2);
        allZsymSel.set('entitydim', 2).set('condition', 'allvertices');
        beamgeom.runCurrent;

        PMLZsymSel = beamgeom.create('PMLZsymSel', 'DifferenceSelection');
        PMLZsymSel.set('entitydim', 2);
        PMLZsymSel.set('add', 'allZsymSel').set('subtract', 'structZsymSel');
        beamgeom.runCurrent;
        P.bndSel.PMLZsym = readSel(model, P.geomname, 'PMLZsymSel');
    end

    PMLxoutSel = beamgeom.create('PMLxoutSel', 'BoxSelection');
    PMLxoutSel.set('xmin', xP1 - d2).set('xmax', xP1 + d2);
    PMLxoutSel.set('ymin', -d2).set('ymax', yP1 + d2);
    PMLxoutSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
    PMLxoutSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;

    PMLyoutSel = beamgeom.create('PMLyoutSel', 'BoxSelection');
    PMLyoutSel.set('xmin', xB1 - d2).set('xmax', xP1 + d2);
    PMLyoutSel.set('ymin', yP1 - d2).set('ymax', yP1 + d2);
    PMLyoutSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
    PMLyoutSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;

    P.bndSel.PMLcurv  = [readSel(model, P.geomname, 'PMLxoutSel'), ...
                         readSel(model, P.geomname, 'PMLyoutSel')];
    P.bndSel.beamXend = P.bndSel.PMLcurv;
else
    P.bndSel.beamXend = [P.bndSel.shieldXend, P.bndSel.shieldYend];
end

%% Audit every emitted selection
% The legacy PML selections resolve to empty and nothing notices (see the KNOWN
% BUG note in the header). Warn loudly here instead.
assertNonEmptySel('P.domSel.beam',   P.domSel.beam);
assertNonEmptySel('P.domSel.shield', P.domSel.shield);
assertNonEmptySel('P.bndSel.beamXsym', P.bndSel.beamXsym);
assertNonEmptySel('P.bndSel.beamYsym', P.bndSel.beamYsym);
assertNonEmptySel('P.bndSel.beamXend', P.bndSel.beamXend);
assertNonEmptySel('P.bndSel.shieldXend', P.bndSel.shieldXend);
assertNonEmptySel('P.bndSel.shieldYend', P.bndSel.shieldYend);
if symZOn
    assertNonEmptySel('P.bndSel.beamZsym', P.bndSel.beamZsym);
end
if usePML
    assertNonEmptySel('P.domSel.PML',      P.domSel.PML);
    assertNonEmptySel('P.bndSel.PMLYsym',  P.bndSel.PMLYsym);
    assertNonEmptySel('P.bndSel.PMLcurv',  P.bndSel.PMLcurv);
    if symZOn
        assertNonEmptySel('P.bndSel.PMLZsym', P.bndSel.PMLZsym);
    end
    if numel(P.domSel.PML) ~= 3
        warning('BuildNanobeamCrossShieldFEM:PMLDomainCount', ...
            ['The PML frame resolved to %d domain(s), not 3. The corner block ' ...
             'must stay separate so COMSOL stretches it in BOTH x and y; if ' ...
             'the blocks merged, the corner absorbs in one direction only.'], ...
            numel(P.domSel.PML));
    end
end

end

% -------------------------------------------------------------------------

function P = readCrossShieldParams(P)
%READCROSSSHIELDPARAMS Fill in defaults and validate, before COMSOL sees anything.

required = {'aShield','hShield','wShield','nShieldX','nShieldY', ...
            'w','th','xsect','geomHalf','beamLenHalf', ...
            'solveMech','mevenx','meveny','mevenz'};
missing = required(~isfield(P, required));
if ~isempty(missing)
    error('BuildNanobeamCrossShieldFEM:missingFields', ...
        'Missing required P field(s): %s.', strjoin(missing, ', '));
end

validateattributes(P.aShield,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.aShield');
validateattributes(P.hShield,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.hShield');
validateattributes(P.wShield,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.wShield');
validateattributes(P.nShieldX, {'numeric'}, {'scalar','real','finite','integer','positive'}, mfilename, 'P.nShieldX');
validateattributes(P.nShieldY, {'numeric'}, {'scalar','real','finite','integer','positive'}, mfilename, 'P.nShieldY');
validateattributes(P.w,        {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.w');
validateattributes(P.th,       {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.th');

% Hard error, exactly as buildCrossStrip.m:386: at h >= a the cross arms reach
% the cell boundary and the solid separates into disconnected corner islands.
% COMSOL builds it happily.
if P.hShield >= P.aShield
    error('BuildNanobeamCrossShieldFEM:armTooLong', ...
        ['P.hShield = %g nm must be strictly less than P.aShield = %g nm. At ' ...
         'h >= a the cross arms reach the cell boundary and the shield ' ...
         'separates into disconnected corner islands.'], ...
        P.hShield*1e9, P.aShield*1e9);
end
if P.wShield >= P.hShield
    warning('BuildNanobeamCrossShieldFEM:armWidthExceedsLength', ...
        ['P.wShield = %g nm is not smaller than P.hShield = %g nm, so the ' ...
         'cross void is degenerate (a square at w = h). The fillet selections ' ...
         'assume a genuine plus shape.'], P.wShield*1e9, P.hShield*1e9);
end

% --- unsupported configurations, v1 -------------------------------------
if ~strcmp(P.xsect,'rect')
    error('BuildNanobeamCrossShieldFEM:xsectNotRect', ...
        ['P.xsect = ''%s'' is not supported. ''tri''/''isoFit'' force ' ...
         'P.mevenz = 0 and make the beam an intersection with a prism; the ' ...
         'junction to a prismatic cross slab is undefined. Use ''rect''.'], ...
        P.xsect);
end
if isfield(P,'asymCav') && P.asymCav
    error('BuildNanobeamCrossShieldFEM:asymCavUnsupported', ...
        ['P.asymCav = 1 is not supported in v1: it needs a mirrored shield ' ...
         'and three more PML blocks at -x. addCrossShield already takes an ' ...
         '(xEdge, sgn) pair so that is a sign flip, not a rewrite.']);
end
if isfield(P,'solveOpt') && P.solveOpt
    error('BuildNanobeamCrossShieldFEM:solveOptUnsupported', ...
        ['P.solveOpt = 1 is not supported in v1: the air cylinder (radius ' ...
         '~2*lambda + w/2) would swallow the shield and the PML frame.']);
end
if ~P.solveMech
    error('BuildNanobeamCrossShieldFEM:solveMechRequired', ...
        'P.solveMech = 1 is required - this builder exists to get a mechanical Q.');
end
if abs(P.meveny) ~= 1
    error('BuildNanobeamCrossShieldFEM:menevyRequired', ...
        ['abs(P.meveny) must be 1. There is no full-y build path in this ' ...
         'pipeline, and with meveny = 0 SetupNanobeamFEM puts the y = 0 face ' ...
         'in neither the symmetric nor the antisymmetric list, leaving a ' ...
         'half beam with a FREE mid-plane - a different structure.']);
end

% --- optional fields ----------------------------------------------------
if ~isfield(P,'shieldPadLen'),     P.shieldPadLen     = 0;       end
if ~isfield(P,'shieldFillet'),     P.shieldFillet     = 'none';  end
if ~isfield(P,'r1Shield'),         P.r1Shield         = 0;       end
if ~isfield(P,'r2Shield'),         P.r2Shield         = 0;       end
if ~isfield(P,'shieldMinFeature'), P.shieldMinFeature = 50e-9;   end
if ~isfield(P,'shieldHmax') || isempty(P.shieldHmax)
    % Three elements across a ligament - the ligaments are the narrowest solid
    % feature in the whole model and they carry the entire shield response. The
    % th/3 term only bites for a shield whose ligaments are wider than the slab
    % is thick, where the thickness becomes the limiting dimension instead.
    P.shieldHmax = min((P.aShield - P.hShield)/3, P.th/3);
end

validateattributes(P.shieldPadLen, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.shieldPadLen');
validateattributes(P.r1Shield,     {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r1Shield');
validateattributes(P.r2Shield,     {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r2Shield');
validateattributes(P.shieldHmax,   {'numeric'}, {'scalar','real','finite','positive'},    mfilename, 'P.shieldHmax');

if ~any(strcmp(P.shieldFillet, {'none','analytic','legacy'}))
    error('BuildNanobeamCrossShieldFEM:badFilletMode', ...
        'P.shieldFillet must be ''none'', ''analytic'' or ''legacy'', not ''%s''.', ...
        P.shieldFillet);
end

if isfield(P,'solveMechPML') && P.solveMechPML
    if ~isfield(P,'PMLLen')
        error('BuildNanobeamCrossShieldFEM:noPMLLen', ...
            'P.solveMechPML = 1 requires P.PMLLen [m].');
    end
    validateattributes(P.PMLLen, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.PMLLen');
    if ~isfield(P,'PMLLenY') || isempty(P.PMLLenY)
        P.PMLLenY = P.PMLLen;
    end
    validateattributes(P.PMLLenY, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.PMLLenY');
end

% The shield must be at least as tall as the half beam, or the beam end face
% hangs over the top of the shield and the union is nonsense.
if P.nShieldY*P.aShield < P.w/2
    error('BuildNanobeamCrossShieldFEM:shieldNarrowerThanBeam', ...
        ['The shield reaches y = %.1f nm but the half beam is %.1f nm wide. ' ...
         'Increase P.nShieldY or P.aShield.'], ...
        P.nShieldY*P.aShield*1e9, P.w/2*1e9);
end
end

% -------------------------------------------------------------------------

function shieldTag = addCrossShield(beamgeom, P, xEdge, sgn, thi)
%ADDCROSSSHIELD Build one 2D cross-cell shield slab and return its 3D tag.
%
%   xEdge  x coordinate of the beam end this shield attaches to [m]
%   sgn    +1 to grow towards +x, -1 to mirror towards -x (asymCav, not wired
%          up in v1 - the argument exists so that adding it is a sign flip)
%
% The slab is one rectangle spanning the clamp pad AND the cells; the pad is
% simply extra length at the attached end, with no cross voids in it, so no
% union is needed to assemble the footprint.

if sgn > 0
    sfx = 'R';
else
    sfx = 'L';
end

as = P.aShield;
hs = P.hShield;
ws = P.wShield;
nx = P.nShieldX;
ny = P.nShieldY;

baseLen = P.shieldPadLen + nx*as;
xS0     = xEdge + sgn*P.shieldPadLen;

wpTag = ['shieldWP', sfx];
shieldWP = beamgeom.feature.create(wpTag, 'WorkPlane');
shieldWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -thi/2);
wpGeom = shieldWP.geom;

baseTag = ['r_shieldBase', sfx];
rBase = wpGeom.feature.create(baseTag, 'Rectangle');
rBase.label('Shield footprint');
rBase.set('type', 'solid').set('base', 'corner');
if sgn > 0
    rBase.set('pos', [xEdge 0]);
else
    rBase.set('pos', [xEdge - baseLen, 0]);
end
rBase.set('size', [baseLen, ny*as]);

xcentres = xS0 + sgn*(((1:nx) - 0.5)*as);
ycentres =           ((1:ny) - 0.5)*as;

subTags = addCrossVoidArray(wpGeom, sfx, xcentres, ycentres, hs, ws);

% strjoin on a cellstr, not strcat in a loop: strcat(char,cell) returns a cell
% and also strips trailing whitespace from char arguments (buildCrossStrip.m
% :212-213). For very large arrays a Difference feature with an input2 cellstr
% avoids the formula string entirely and is equivalent.
composeTag = ['shieldCo', sfx];
composeFeat = wpGeom.feature.create(composeTag, 'Compose');
composeFeat.selection('input').set([{baseTag}, subTags]);
composeFeat.set('formula', [baseTag, '-', strjoin(subTags, '-')]);
composeFeat.set('intbnd', false);

% One runCurrent for the whole array, not one per cell. buildCrossUnitCell runs
% the sequence after almost every feature, which for nx*ny cells would mean
% O(nx*ny) full geometry rebuilds (buildCrossStrip.m:189-192).
beamgeom.runCurrent;

% Fillets, if any. Default 'none' emits nothing at all - COMSOL's Fillet errors
% on radius 0 anyway.
if ~strcmp(P.shieldFillet,'none') && (P.r1Shield > 0 || P.r2Shield > 0)
    selWidth = 10e-9;
    for i = 1:nx
        for j = 1:ny
            addCrossFillets(wpGeom, sfx, i, j, xcentres(i), ycentres(j), ...
                            hs, ws, P.r1Shield, P.r2Shield, ...
                            P.shieldFillet, selWidth);
        end
    end
    beamgeom.runCurrent;
end

shieldTag = ['shieldExt', sfx];
shieldExt = beamgeom.feature.create(shieldTag, 'Extrude');
shieldExt.set('distance', thi);
shieldExt.selection('input').set({wpTag});
beamgeom.runCurrent;
end

% -------------------------------------------------------------------------

function subTags = addCrossVoidArray(wpGeom, sfx, xcentres, ycentres, hs, ws)
%ADDCROSSVOIDARRAY Two crossed rectangles per cell over a 2D (i,j) grid.
%
% Ported from buildCrossStrip.m:194-210, with a second index. Every feature is
% created first and the sequence is run once by the caller.

nx = numel(xcentres);
ny = numel(ycentres);

subTags = cell(1, 2*nx*ny);
k = 0;
for i = 1:nx
    for j = 1:ny
        tagH = sprintf('cs%s_%d_%d_h', sfx, i, j);
        tagV = sprintf('cs%s_%d_%d_v', sfx, i, j);

        recH = wpGeom.feature.create(tagH, 'Rectangle');
        recH.set('base', 'center');
        recH.set('pos',  [xcentres(i) ycentres(j)]);
        recH.set('size', [hs ws]);

        recV = wpGeom.feature.create(tagV, 'Rectangle');
        recV.set('base', 'center');
        recV.set('pos',  [xcentres(i) ycentres(j)]);
        recV.set('size', [ws hs]);

        subTags{k+1} = tagH;
        subTags{k+2} = tagV;
        k = k + 2;
    end
end
end

% -------------------------------------------------------------------------

function addCrossFillets(wpGeom, sfx, i, j, xc, yc, hs, ws, r1, r2, mode, selWidth)
%ADDCROSSFILLETS Per-cell DiskSelection + Fillet pairs.
%
%   mode = 'legacy'    the annuli of buildCrossUnitCell.m:129-148, unchanged.
%                      disksel1 = [w/2, 1.5*w], disksel2 = [h/2 +/- 5nm].
%                      For the legacy cross parameters disksel2 selects nothing
%                      (the tips are at sqrt(h^2+w^2)/2, not h/2) and disksel1
%                      catches all 12 own-cell corners plus 8 foreign ones.
%   mode = 'analytic'  thin annuli on the actual vertex radii, so r1 rounds the
%                      four inner concave corners and r2 the eight arm tips.
%
% A zero radius is skipped: COMSOL's Fillet errors on radius 0.

switch mode
    case 'legacy'
        rin1 = ws/2;              rout1 = ws*(3/2);
        rin2 = hs/2 - selWidth/2; rout2 = hs/2 + selWidth/2;
    case 'analytic'
        Rin  = ws/sqrt(2);                % 4 inner concave corners
        Rtip = sqrt(hs^2 + ws^2)/2;       % 8 arm-tip corners
        rin1 = Rin  - selWidth/2;  rout1 = Rin  + selWidth/2;
        rin2 = Rtip - selWidth/2;  rout2 = Rtip + selWidth/2;
    otherwise
        error('BuildNanobeamCrossShieldFEM:badFilletMode', ...
            'Unknown fillet mode ''%s''.', mode);
end

if r1 > 0
    selTag = sprintf('cs%s_%d_%d_ds1', sfx, i, j);
    ds1 = wpGeom.feature.create(selTag, 'DiskSelection');
    ds1.set('entitydim', 0);
    ds1.set('posx', xc).set('posy', yc);
    ds1.set('r', rout1).set('rin', rin1);
    ds1.set('condition', 'allvertices');

    fil1 = wpGeom.feature.create(sprintf('cs%s_%d_%d_fil1', sfx, i, j), 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(selTag);
end

if r2 > 0
    selTag = sprintf('cs%s_%d_%d_ds2', sfx, i, j);
    ds2 = wpGeom.feature.create(selTag, 'DiskSelection');
    ds2.set('entitydim', 0);
    ds2.set('posx', xc).set('posy', yc);
    ds2.set('r', rout2).set('rin', rin2);
    ds2.set('condition', 'allvertices');

    fil2 = wpGeom.feature.create(sprintf('cs%s_%d_%d_fil2', sfx, i, j), 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(selTag);
end
end

% -------------------------------------------------------------------------

function auditCrossShield(P, wid, thi, usePML)
%AUDITCROSSSHIELD Pure-MATLAB, warn-only checks run before anything is built.
%
% Nothing here changes the geometry. The point is that every compromise in this
% design is numeric and knowable in advance, so none of them should be silent.

as = P.aShield;
hs = P.hShield;
ws = P.wShield;

wall     = as - hs;         % solid between the voids of adjacent cells
nodeHalf = (as - ws)/2;     % half-width of a node, measured from a cell corner
finThk   = wall/2;          % half wall left at a free cut through a ligament

% --- lithography --------------------------------------------------------
if exist('isCrossFabricable','file') == 2
    [ok, margin, why] = isCrossFabricable(as, hs, ws, thi, ...
                                          P.shieldMinFeature, ...
                                          P.r1Shield, P.r2Shield);
    if ~ok
        warning('BuildNanobeamCrossShieldFEM:notFabricable', ...
            ['Cross cell fails isCrossFabricable at minFeature = %g nm: ' ...
             'tightest rule is ''%s'', short by %.1f nm. The binding ' ...
             'constraint is normally the inter-cell wall a-h = %.1f nm.'], ...
            P.shieldMinFeature*1e9, why{1}, -margin(1)*1e9, wall*1e9);
    end
else
    warning('BuildNanobeamCrossShieldFEM:noFabCheck', ...
        ['isCrossFabricable.m is not on the path (it lives in ' ...
         'omc-comsol-chang/), so the lithography audit was skipped.']);
end

% --- how the beam actually lands on the shield --------------------------
if P.shieldPadLen <= 0
    if nodeHalf < wid/2
        fracOnFin = (wid/2 - nodeHalf)/(wid/2);
        warning('BuildNanobeamCrossShieldFEM:beamOverhangsNode', ...
            ['P.shieldPadLen = 0 and the node is narrower than the beam: ' ...
             'node half-width (a-w)/2 = %.1f nm vs beam half-width %.1f nm. ' ...
             '%.0f%% of the beam end half-face lands on the %.1f nm ' ...
             'half-wall fin instead of on solid material. Set ' ...
             'P.shieldPadLen > 0 (about a/2) for a proper clamp, or choose ' ...
             'a - w >= %.1f nm.'], ...
            nodeHalf*1e9, wid/2*1e9, 100*fracOnFin, finThk*1e9, wid*1e9);
    end
    warning('BuildNanobeamCrossShieldFEM:freeEdgeFins', ...
        ['P.shieldPadLen = 0 leaves the shield''s x = xS0 edge cut through ' ...
         'the y-running ligaments, so %d free-standing fins of thickness ' ...
         '(a-h)/2 = %.1f nm hang off it. They are supported on two edges and ' ...
         'will contribute low-frequency flapping modes near the free edge. ' ...
         'P.shieldPadLen > 0 bonds them to the clamp slab.'], ...
        P.nShieldY, finThk*1e9);
end

% --- fillet selection bleed --------------------------------------------
% The cross void is 4-fold symmetric, so all four orthogonal neighbours
% contribute the same two nearest vertices.
if ~strcmp(P.shieldFillet,'none') && (P.r1Shield > 0 || P.r2Shield > 0)
    selWidth = 10e-9;
    Rin      = ws/sqrt(2);
    Rtip     = sqrt(hs^2 + ws^2)/2;
    % Nearest vertex belonging to an orthogonal neighbour one pitch away. The
    % diagonal neighbours sit at sqrt(2)*a and never reach.
    Rforeign = sqrt((as - hs/2)^2 + (ws/2)^2);

    switch P.shieldFillet
        case 'legacy'
            rout1 = ws*(3/2);
            if Rforeign <= rout1 && P.r1Shield > 0
                warning('BuildNanobeamCrossShieldFEM:filletSelectionBleed', ...
                    ['Legacy disksel1 spans [%.1f, %.1f] nm and the nearest ' ...
                     'FOREIGN vertex sits at %.1f nm, so each cell''s r1 ' ...
                     'fillet rounds 8 corners belonging to its four ' ...
                     'orthogonal neighbours - twice the bleed ' ...
                     'buildCrossStrip.m documents for a 1D strip. With %d ' ...
                     'cells that is %d Fillet features rounding ' ...
                     'already-rounded corners. EXPECT COMSOL TO REJECT THIS. ' ...
                     'Use P.shieldFillet = ''analytic'' or ''none''.'], ...
                    ws/2*1e9, rout1*1e9, Rforeign*1e9, ...
                    P.nShieldX*P.nShieldY, 2*P.nShieldX*P.nShieldY);
            end
            if P.r2Shield > 0 && ...
               (Rtip < hs/2 - selWidth/2 || Rtip > hs/2 + selWidth/2)
                warning('BuildNanobeamCrossShieldFEM:filletR2NoOp', ...
                    ['Legacy disksel2 spans [%.1f, %.1f] nm but the arm-tip ' ...
                     'corners sit at %.1f nm, so it selects nothing and ' ...
                     'P.r2Shield = %g nm has no effect. This reproduces ' ...
                     'buildCrossUnitCell exactly; it is reported, not ' ...
                     'corrected.'], ...
                    (hs/2-selWidth/2)*1e9, (hs/2+selWidth/2)*1e9, ...
                    Rtip*1e9, P.r2Shield*1e9);
            end
        case 'analytic'
            marginTip = Rforeign - (Rtip + selWidth/2);
            marginIn  = Rforeign - (Rin  + selWidth/2);
            if marginTip <= 0 || marginIn <= 0
                warning('BuildNanobeamCrossShieldFEM:analyticFilletBleed', ...
                    ['The analytic fillet annuli reach a neighbouring cell: ' ...
                     'nearest foreign vertex at %.1f nm vs tip annulus outer ' ...
                     'radius %.1f nm (margin %.1f nm) and inner-corner ' ...
                     'annulus outer radius %.1f nm (margin %.1f nm). Reduce ' ...
                     'P.wShield relative to P.aShield.'], ...
                    Rforeign*1e9, (Rtip+selWidth/2)*1e9, marginTip*1e9, ...
                    (Rin+selWidth/2)*1e9, marginIn*1e9);
            end
    end

    % Fillet radii must fit the edges they consume (isCrossFabricable's rules).
    if P.r2Shield > ws/2
        warning('BuildNanobeamCrossShieldFEM:filletTooLarge', ...
            'P.r2Shield = %.1f nm exceeds w/2 = %.1f nm; the arm end cannot take it.', ...
            P.r2Shield*1e9, ws/2*1e9);
    end
    if P.r1Shield + P.r2Shield > (hs - ws)/2
        warning('BuildNanobeamCrossShieldFEM:filletTooLarge', ...
            'r1+r2 = %.1f nm exceeds the arm side (h-w)/2 = %.1f nm.', ...
            (P.r1Shield+P.r2Shield)*1e9, (hs-ws)/2*1e9);
    end
end

% --- PML thickness against the mechanical wavelength --------------------
if usePML && isfield(P,'freq') && P.freq > 0
    v = [];
    if isfield(P,'D') && numel(P.D) >= 10 && isfield(P,'rho')
        v = sqrt(P.D(10)/P.rho);                 % c44 -> shear
    elseif isfield(P,'E') && isfield(P,'nu') && isfield(P,'rho')
        v = sqrt(P.E/(2*(1+P.nu))/P.rho);
    end
    if ~isempty(v)
        lam = v/P.freq;
        if P.PMLLen < 0.5*lam || P.PMLLenY < 0.5*lam
            warning('BuildNanobeamCrossShieldFEM:PMLThin', ...
                ['PML is thin: %.2f um (+x) / %.2f um (+y) against a shear ' ...
                 'wavelength of %.2f um at %.2f GHz. Aim for about one ' ...
                 'wavelength and confirm with a Q vs PMLLen sweep.'], ...
                P.PMLLen*1e6, P.PMLLenY*1e6, lam*1e6, P.freq/1e9);
        end
        if P.PMLLen > 4*lam
            warning('BuildNanobeamCrossShieldFEM:PMLThick', ...
                ['PML is %.1f wavelengths thick (%.2f um vs %.2f um). The ' ...
                 'frame plan area is %.1f um^2 against a shield plan area of ' ...
                 '%.1f um^2, so most of the mesh will be spent on the ' ...
                 'absorber.'], ...
                P.PMLLen/lam, P.PMLLen*1e6, lam*1e6, ...
                (P.PMLLen*(P.nShieldY*as + P.PMLLenY) + ...
                 (P.nShieldX*as + P.shieldPadLen)*P.PMLLenY)*1e12, ...
                (P.nShieldX*as)*(P.nShieldY*as)*1e12);
        end
    end
end

% --- mesh cost ----------------------------------------------------------
if P.shieldHmax > wall/2
    warning('BuildNanobeamCrossShieldFEM:shieldMeshCoarse', ...
        ['P.shieldHmax = %.1f nm puts fewer than two elements across the ' ...
         '%.1f nm ligament. The ligaments carry the whole shield response; ' ...
         'under-resolving them makes the Q meaningless.'], ...
        P.shieldHmax*1e9, wall*1e9);
end
nElemEst = P.shield.solidFrac * (P.nShieldX*as) * (P.nShieldY*as) * ...
           (thi/2) / (P.shieldHmax^3/6);
if nElemEst > 2e6
    warning('BuildNanobeamCrossShieldFEM:shieldMeshCost', ...
        ['Rough tet-count estimate for the shield alone is %.1e elements at ' ...
         'hmax = %.1f nm (solid fraction %.2f). P.max_dof is %g, so ' ...
         'SolveNanobeamFEM''s mAdjMesh loop will coarsen until the ligaments ' ...
         'are unresolved. Reduce P.nShieldX/P.nShieldY, or widen a-h.'], ...
        nElemEst, P.shieldHmax*1e9, P.shield.solidFrac, ...
        localGetField(P,'max_dof',NaN));
end
end

% -------------------------------------------------------------------------

function v = localGetField(S, name, dflt)
%LOCALGETFIELD Field value or a default, for warning text only.
if isfield(S, name)
    v = S.(name);
else
    v = dflt;
end
end

% -------------------------------------------------------------------------

function inds = readSel(model, geomname, tag)
%READSEL Resolve a geometry selection feature to a row vector of entity indices.
inds = double(model.selection([geomname,'_',tag]).inputEntities());
inds = reshape(inds, 1, []);
end

% -------------------------------------------------------------------------

function assertNonEmptySel(name, inds)
%ASSERTNONEMPTYSEL Warn when a selection this builder promises resolved to nothing.
%
% This exists because of BuildNanobeamFEM.m:508/518/538, where a box that is
% too small silently yields an empty selection, a boundary condition is applied
% to nothing, and the resulting Q looks perfectly plausible.
if isempty(inds)
    warning('BuildNanobeamCrossShieldFEM:emptySelection', ...
        ['%s resolved to NO entities. Any boundary condition, material or ' ...
         'mesh rule hung on it will be applied to nothing. Check the ' ...
         'geometry in the GUI before trusting anything downstream.'], name);
end
end
