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
%   OPTICAL P FIELDS  (read only when P.solveOpt = 1)
%       P.airrad       radius of the air half-cylinder wrapped round the beam
%                      [m]. REQUIRED. The house convention, from
%                      test_nanobeamRectFEM.m:105 and five sibling scripts, is
%                      2*P.lambda + P.w/2; for lambda = 1550 nm and w = 750 nm
%                      that is 3475 nm.
%       P.airCylLen    +x extent of that cylinder [m].
%                      Default min(xB1, xS0 - P.lambda). See WHY THE AIR
%                      CYLINDER STOPS SHORT below. Hard error if the cylinder
%                      envelope touches the shield or the PML.
%       P.lambda       target vacuum wavelength [m]. REQUIRED (also used for
%                      the default P.airCylLen).
%       P.nbeam        beam refractive index. REQUIRED by SetupNanobeamFEM.
%       P.oevenx/P.oeveny/P.oevenz   optical symmetry, +/-1 (0 for oevenx only,
%                      which puts a scattering BC on x = 0 instead). P.oevenz
%                      must agree with P.mevenz on whether z = 0 is a mirror
%                      plane at all - see the hard error in
%                      readCrossShieldParams.
%
%   WHY THE AIR CYLINDER STOPS SHORT  (the reason P.solveOpt used to be barred)
%     BuildNanobeamFEM.m:202-230 revolves an air half-cylinder of radius
%     P.airrad over the WHOLE beam length. Ported verbatim that cylinder does
%     genuinely swallow this geometry: it reaches |z| = P.airrad = 3475 nm
%     against a slab that is only +/-250 nm thick, and y = 3475 nm against a
%     shield that starts at x = xS0 and runs to y = yS1 = 8000 nm. A cylinder
%     spanning x in [0, xS1] would enclose most of the shield in air, and the
%     shield would stop being a shield.
%
%     The cavity mode, however, lives at x = 0 and is exponentially attenuated
%     by the hole mirror long before the shield. So the cylinder is scoped to
%     the CAVITY, not to the structure: x in [0, P.airCylLen] with
%     P.airCylLen < xS0, and the shield and the PML sit in bare vacuum that is
%     simply not meshed. Optically the beam beyond P.airCylLen is clad in the
%     default PEC of the EM interface, which is harmless exactly because the
%     mirror has already killed the field there - stage S4b in
%     test_nanobeamCrossShieldFEM.m is the convergence check that proves it.
%
%     The cylinder does NOT end flush at xS0. A flush end cap would share the
%     x = xS0 plane with the shield's front face, so the part of the cap inside
%     |y| < airrad, |z| < th/2 would become an INTERIOR air|diamond boundary
%     and land in P.bndSel.cylXend, which SetupNanobeamFEM.m:390 feeds to a
%     Scattering BC. A scattering BC on an interior face is an absorber buried
%     in the model. Hence the default leaves one free-space wavelength of
%     vacuum between the cap and the shield.
%
%     Worked numbers for test_nanobeamCrossShieldFEM.m (all nm):
%         xB1 = beamLenHalf = 8951.8,  shieldPadLen = 0  ->  xS0 = 8951.8
%         xS1 = 8951.8 + 6*1600 = 18551.8,  yS1 = 5*1600 = 8000
%         xP1 = 20551.8,  yP1 = 10000
%         airrad = 2*1550 + 375 = 3475
%         airCylLen = min(8951.8, 8951.8 - 1550) = 7401.8
%       cylinder envelope   x [0, 7401.8]   y [0, 3475]   z [0, 3475]
%       shield envelope     x [8951.8, 18551.8]  y [0, 8000]  z [0, 250]
%         -> separated in x by 1550.0 nm.  y and z OVERLAP, so x is the only
%            separation and it is the thing the hard error guards.
%       PML +x arm          x [18551.8, 20551.8] -> separated in x by 11150.0
%       PML +y arm          x [8951.8, 18551.8]  y [8000, 10000]
%         -> separated in x by 1550.0 nm AND in y by 4525.0 nm
%       PML corner          separated in x by 11150.0 and in y by 4525.0
%       The cap at x = 7401.8 falls in solid beam between holes 15 (ends
%       7198.8) and 16 (starts 7530.8), so 15 of the 18 half-beam holes - the
%       6 defect holes and 9 mirror periods - are inside the air.
%
%   THE END CAP SPLITS THE BEAM, ON PURPOSE
%     Compose keeps interior boundaries (it must: the beam and the air have to
%     stay separate domains so they can carry different materials), so the end
%     cap is imprinted on the beam as well and the beam becomes TWO domains,
%     x < airCylLen and x > airCylLen. Both are picked up by the 'inside' boxes
%     that build P.domSel.beam, so nothing downstream notices. What does notice
%     is P.bndSel.cylXend: the imprinted beam cross-section lies in the same
%     x = airCylLen plane as the cap, so it is removed with a
%     DifferenceSelection, leaving only faces with vacuum on the outside.
%
%     One second-order consequence to be aware of rather than to fix here.
%     CalcGOM.m:78-84 builds its moving-boundary integration surface as an
%     AdjacentSelection on the geometry feature 'beamSel' - which this builder
%     provides - and then removes only the three symmetry planes. That surface
%     therefore includes the cap imprint at x = airCylLen and the beam|shield
%     interface at x = xS0, both of which are diamond|diamond and should
%     contribute nothing to a moving-boundary term. They do not, in practice,
%     because both sit 9+ mirror periods from the defect where |E|^2 is already
%     down by orders of magnitude - the same argument that justifies truncating
%     the air in the first place. Stage S4b in the test script is what turns
%     that argument into a number.
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
%     and, only when P.solveOpt = 1, the optical half of the contract that
%     SetupNanobeamFEM.m:176 and :368-419 consume:
%       P.domSel.cyl      the air domain(s)  (air material, emw selection)
%       P.bndSel.cylXsym  x = 0 plane, beam AND air faces - the full optical
%                         mirror plane, not just the solid part
%       P.bndSel.cylXend  x = airCylLen air end cap, EXTERIOR faces only (the
%                         imprinted beam cross-section is differenced out)
%       P.bndSel.cylYsym  y = 0 plane, beam + pad + shield + air faces
%       P.bndSel.cylZsym  z = 0 plane, beam + pad + shield + air faces
%                         (if symZOn)
%       P.bndSel.cylCurv  the curved outer surface of the air cylinder
%       P.airCyl          struct of computed cylinder extents
%
%     Note that P.bndSel.beamZsym and P.bndSel.cylZsym DIFFER once the air
%     exists: the mechanical symmetry BC must not be hung on an air face (no
%     Solid Mechanics is solved there), while the optical PMC/PEC must cover
%     the whole plane. beamZsym is cylZsym with the air faces removed.
%
%   NOT SUPPORTED IN v1 (all hard errors, see the checks below)
%       P.asymCav      - needs a mirrored shield and three more PML blocks at
%                        -x. addCrossShield takes an (xEdge, sgn) pair so the
%                        left end is a sign flip rather than a rewrite. The
%                        air cylinder would also need mirroring to -x.
%       P.xsect ~= 'rect' - 'tri'/'isoFit' force P.mevenz = 0 and make the beam
%                        an intersection with a prism; the junction to a
%                        prismatic cross slab is undefined.
%       P.meveny == 0  - there is no full-y build path in this pipeline.
%       abs(P.oevenz) ~= abs(P.mevenz) when P.solveOpt = 1 - there is a single
%                        z = 0 cut shared by both physics, so either both want
%                        it or neither does.
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
%   THE ONE MESH NODE, AND WHY THE COMBINED RUN IS EXPENSIVE
%     This is a geometry builder, so it cannot fix it - but a caller turning
%     P.solveOpt on needs to know. SolveNanobeamFEM creates a SINGLE mesh node
%     ('mesh') and both studies share it. The optical branch runs first with
%     custom hmax = P.lambda/5 applied GLOBALLY (SolveNanobeamFEM.m:131) and no
%     domain overrides, then the mechanical branch turns that custom size OFF
%     (:307) and re-runs the mesh with hmax = P.shieldHmax on P.domSel.shield.
%
%     For test_nanobeamCrossShieldFEM.m's parameters the two meshes are not
%     remotely the same object in the shield: 310 nm gives about 1.7e3 tets
%     there, 33.3 nm gives about 1.4e6. The optical study therefore solves on a
%     mesh in which the shield is barely discretised (harmless, the field there
%     is zero) and the mechanical study solves on a mesh in which the air
%     cylinder has no size rule at all. CalcGOM.m:124-155 then Joins the two
%     datasets, which requires a common mesh.
%
%     Per-study DOF counts stay separate, because each physics is selected on
%     its own domains: emw gets beam + pad + shield + air, smech gets beam +
%     pad + shield + PML, and neither sees the other's volume. So the optical
%     study is cheap (order 5e5 DOF) and adding it does NOT raise the
%     mechanical DOF count. The mechanical side is the binding constraint, and
%     it is already binding without optics - see the estimate in
%     auditCrossShield and the S3 note in the test script.
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

% Air half-cylinder around the CAVITY only. See WHY THE AIR CYLINDER STOPS
% SHORT in the header. useOpt is the single gate for every optical addition in
% this file, so P.solveOpt = 0 runs exactly the code it ran before.
useOpt = isfield(P,'solveOpt') && P.solveOpt;
if useOpt
    airrad = P.airrad;
    xC1    = P.airCylLen;           % +x extent of the cylinder
else
    airrad = 0;
    xC1    = 0;
end

% z symmetry is only meaningful for the rectangular cross-section, and v1 is
% rect only, so the guard collapses to the symmetry flags. With optics on BOTH
% physics have to agree that z = 0 is a mirror plane, because there is only one
% geometric cut; readCrossShieldParams hard-errors if they disagree, so the
% && below can never be a silent downgrade.
symZOn = abs(P.mevenz) > 0 && (~useOpt || abs(P.oevenz) > 0);
if symZOn
    zLoEff = 0;                     % solid slab reaches z = 0 from below
    zLoCyl = 0;                     % air cylinder likewise
else
    zLoEff = -thi/2;
    zLoCyl = -airrad;
end

% The z-symmetry cut block has to span the FULL z reach of the model, not just
% the slab. BuildNanobeamFEM.m:293-295 has this as commented-out intent
% (symZthList(end+1) = P.airrad); leaving it at thi/2 would cut the slab and
% leave the z < 0 half of the air cylinder standing, so the z = 0 optical
% mirror plane would be a plane through the middle of a solid air region
% instead of an exterior boundary - no faces at z = 0 for the air, and a
% geometry that is not the eighth model the symmetry factors in CalcGOM assume.
zCutLo = max(thi/2, airrad);

if useOpt
    % Hard error, not a warning: the guard this replaces existed precisely to
    % stop the cylinder eating the shield, and a later parameter change must
    % not be able to recreate that silently.
    assertCylClears(P, xC1, airrad, zLoCyl, xB1, xS0, xS1, yS1, xP1, yP1, ...
                    zLoEff, thi, usePML);

    P.airCyl = struct( ...
        'rad',      airrad, ...
        'len',      xC1, ...
        'gapToxS0', xS0 - xC1, ...
        'zLo',      zLoCyl, ...
        'nHolesIn', sum(xpos + hx/2 <= xC1));
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
auditCrossShield(P, wid, thi, usePML, useOpt);

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

%% Air half-cylinder around the cavity
% Ported from BuildNanobeamFEM.m:202-230: an xz work plane carrying a
% rectangle at pos [0 0] of size [len airrad], revolved about axis [1 0] from
% angle1 = -180 to angle2 = 0. That is a HALF cylinder covering y >= 0, which
% is the half the rest of the model lives in.
%
% Two departures from the reference:
%   * the rectangle length is P.airCylLen, not the full beam length, so the
%     cylinder never reaches the shield (see the header);
%   * the input selection of the Revolve is set explicitly. The reference omits
%     it and relies on the default picking the only work plane, and this file
%     already has three.
displayCylStr = '';
if useOpt
    cylWP = beamgeom.feature.create('cyl_cut_wp', 'WorkPlane');
    cylWP.set('planetype', 'quick').set('quickplane', 'xz');
    cylRect = cylWP.geom.feature.create('cyl_cut', 'Rectangle');
    cylRect.set('type', 'solid').set('base', 'corner');
    cylRect.set('pos', [0 0]).set('size', [xC1, airrad]);
    beamgeom.runCurrent;

    airCyl = beamgeom.feature.create('air_cyl', 'Revolve');
    airCyl.selection('input').set({'cyl_cut_wp'});
    airCyl.set('angle1', '-180');
    airCyl.set('angle2', '0');
    airCyl.set('axis', [1,0]).set('pos', [0,0]);
    beamgeom.runCurrent;

    % intbnd is left at the Compose default (interior boundaries KEPT), which
    % is what the reference relies on and what the contract needs: beam and air
    % must stay separate domains or they cannot carry different materials.
    % The price is the end-cap imprint that splits the beam in two - see THE
    % END CAP SPLITS THE BEAM, ON PURPOSE in the header.
    beamAir = beamgeom.feature.create('beamAir', 'Compose');
    beamAir.selection('input').set(finBeamTag);
    beamAir.selection('input').set('air_cyl');
    beamAir.set('formula', [finBeamTag, ' + air_cyl']);
    beamgeom.runCurrent;

    finBeamTag = 'beamAir';
    displayCylStr = sprintf(', air cylinder r = %.0f nm to x = %.0f nm', ...
                            airrad*1e9, xC1*1e9);
end

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
%
% The cut block spans z in [-zCutLo, 0] where zCutLo = max(thi/2, airrad), NOT
% thi/2. With the air cylinder present thi/2 would cut only the slab and leave
% the whole z < 0 half of the cylinder standing (see the zCutLo comment above);
% with P.solveOpt = 0, airrad = 0 and zCutLo == thi/2, so this is the byte-for-
% byte previous behaviour.
%
% yMaxAll = yP1 = 10000 nm and xMaxAll = xP1 cover the cylinder's y and x reach
% (3475 nm and P.airCylLen) with room to spare, so the footprint needs no
% widening for the cylinder - only the z extent does.
if symZOn
    symZWP = beamgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -zCutLo);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [xMinAll 0]).set('size', [xMaxAll - xMinAll, yMaxAll]);
    beamgeom.runCurrent;

    symZPlaneExt = beamgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', zCutLo);
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
      num2str(2*ny), ' physical rows)', displayCylStr, displayPMLStr]);

%% Domain selections
% WHY THE AIR DOMAIN DOES NOT LEAK INTO THE SOLID SELECTIONS
% beamSel, shieldSel, structSel and allSel are all condition 'inside' with
% zmax = thi/2 + d = 260 nm. The air cylinder reaches z = airrad = 3475 nm, so
% it is never 'inside' any of them and needs no explicit exclusion. That is a
% property of the z extent, not of the x/y extents, so do not widen zmax on
% these four boxes. The PML 'intersects' probes sit at x >= xB1, well clear of
% the cylinder's x in [0, airCylLen], so they are unaffected too.
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

% Air domain(s) - BuildNanobeamFEM.m:342-360, with the x extent of the box
% taken from the cylinder rather than from the beam.
%
% The construction is 'inside' box over the cylinder envelope MINUS beamSel.
% beamSel covers BOTH halves of the end-cap-split beam (it is an 'inside' box
% out to xB1 + d), so the subtraction removes the beam whichever way COMSOL
% numbers the two pieces. The beam half at x > airCylLen is not inside the
% cylinder box at all, so it never enters the sum in the first place.
%
% The air is a single connected domain even though the holes look isolated: the
% elliptical holes have semi-axis max(hy)/2 = 289 nm against a half width of
% w/2 = 375 nm, so they do NOT break out through the beam side wall, but they
% are through-holes in z and open onto the z = +thi/2 face, which is inside the
% cylinder. Do not assume numel(P.domSel.cyl) == 1 though - the audit only
% checks that it is non-empty.
if useOpt
    beamCylSel = beamgeom.create('beamCylSel', 'BoxSelection');
    beamCylSel.set('xmin', -d).set('xmax', xC1 + d);
    beamCylSel.set('ymin', -d).set('ymax', airrad + d);
    beamCylSel.set('zmin', zLoCyl - d).set('zmax', airrad + d);
    beamCylSel.set('entitydim', 3).set('condition', 'inside');
    beamgeom.runCurrent;

    cylSel = beamgeom.create('cylSel', 'DifferenceSelection');
    cylSel.set('entitydim', 3).set('add', 'beamCylSel').set('subtract', 'beamSel');
    beamgeom.runCurrent;
    P.domSel.cyl = readSel(model, P.geomname, 'cylSel');
    disp(['Air domain indices: ', num2str(P.domSel.cyl)]);
end

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
%
% zmax = thi/2 + d2 keeps the AIR's y = 0 face out of this selection: the air
% reaches z = airrad, so its y = 0 face has vertices far outside the box and
% 'allvertices' rejects it. That matters - this selection feeds a Solid
% Mechanics symmetry/antisymmetry BC, and there is no Solid Mechanics in the
% air. The optical equivalent, which DOES want the air, is cylYsymSel below.
structYsymSel = beamgeom.create('structYsymSel', 'BoxSelection');
structYsymSel.set('xmin', -d2).set('xmax', xS1 + d2);
structYsymSel.set('ymin', -d2).set('ymax', d2);
structYsymSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
structYsymSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
P.bndSel.beamYsym = readSel(model, P.geomname, 'structYsymSel');

if symZOn
    % z = 0 is the one plane where the air cannot be excluded by a z bound: all
    % faces in it are planar at z = 0, so every vertex of every candidate face
    % is at z = 0 regardless of which domain it belongs to. This box therefore
    % collects the z = 0 faces of beam + pad + shield AND of the air, which is
    % exactly the optical mirror plane (P.bndSel.cylZsym) but NOT the
    % mechanical one.
    structZsymSel = beamgeom.create('structZsymSel', 'BoxSelection');
    structZsymSel.set('xmin', -d2).set('xmax', xS1 + d2);
    structZsymSel.set('ymin', -d2).set('ymax', yS1 + d2);
    structZsymSel.set('zmin', -d2).set('zmax', d2);
    structZsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;

    if ~useOpt
        P.bndSel.beamZsym = readSel(model, P.geomname, 'structZsymSel');
    else
        % Remove the air's z = 0 faces, discriminating on y. There are exactly
        % two families of them and neither overlaps the solid faces in y:
        %
        %   outer air  y in [wid/2, airrad]  - the quarter-disc minus the beam
        %              cross-section, one face spanning x in [0, airCylLen].
        %              Its inner edge sits on the beam side wall at y = wid/2,
        %              its outer edge on the cylinder at y = airrad.
        %   hole air   y in [0, hy_k/2] for each enclosed hole k. hy is at most
        %              max(hy) = 578 nm, so hy/2 <= 289 nm against wid/2 = 375
        %              nm - the 86 nm of headroom that makes this work. It is
        %              checked in readCrossShieldParams, which hard-errors if
        %              the holes are wide enough to reach the side wall.
        %
        % The SOLID z = 0 faces (beam either side of the end-cap split, pad,
        % shield) all run from y = 0 to y = wid/2 or y = yS1, so they fail both
        % boxes: the outer box because of their vertices at y = 0, the hole box
        % because of their vertices at y = wid/2 or beyond.
        %
        % Two chained DifferenceSelections rather than one UnionSelection plus
        % one difference, to stay with the single-string add/subtract form used
        % everywhere else in this file.
        airZ0OutSel = beamgeom.create('airZ0OutSel', 'BoxSelection');
        airZ0OutSel.set('xmin', -d2).set('xmax', xC1 + d2);
        airZ0OutSel.set('ymin', wid/2 - d2).set('ymax', airrad + d2);
        airZ0OutSel.set('zmin', -d2).set('zmax', d2);
        airZ0OutSel.set('entitydim', 2).set('condition', 'allvertices');
        beamgeom.runCurrent;

        airZ0HoleSel = beamgeom.create('airZ0HoleSel', 'BoxSelection');
        airZ0HoleSel.set('xmin', -d2).set('xmax', xC1 + d2);
        airZ0HoleSel.set('ymin', -d2).set('ymax', max(hy)/2 + d2);
        airZ0HoleSel.set('zmin', -d2).set('zmax', d2);
        airZ0HoleSel.set('entitydim', 2).set('condition', 'allvertices');
        beamgeom.runCurrent;

        structZsymNoOut = beamgeom.create('structZsymNoOut', 'DifferenceSelection');
        structZsymNoOut.set('entitydim', 2);
        structZsymNoOut.set('add', 'structZsymSel').set('subtract', 'airZ0OutSel');
        beamgeom.runCurrent;

        structZsymSolid = beamgeom.create('structZsymSolid', 'DifferenceSelection');
        structZsymSolid.set('entitydim', 2);
        structZsymSolid.set('add', 'structZsymNoOut').set('subtract', 'airZ0HoleSel');
        beamgeom.runCurrent;
        P.bndSel.beamZsym = readSel(model, P.geomname, 'structZsymSolid');
    end
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

%% Optical boundary selections
% The exact contract SetupNanobeamFEM consumes:
%   :368/:410  cylXsym  -> PMC (oevenx = +1) / PEC (-1) / scattering (0)
%   :372/:414  cylYsym  -> PMC (oeveny = +1) / PEC (-1)
%   :377/:419  cylZsym  -> PMC (oevenz = +1) / PEC (-1)
%   :390       cylXend and cylCurv -> Scattering BC (always, both required)
%
% The three cyl*sym selections span the WHOLE symmetry plane of the EM region,
% solid and air together, which is why they are separate from the beam*sym
% selections used by Solid Mechanics rather than reusing them. The EM region is
% beam + pad + shield + air, because SetupNanobeamFEM.m:163 takes the optical
% dielectric selection from P.domSel.beam - the same selection Solid Mechanics
% uses - and that one has to contain the shield.
%
% For the PEC branch SetupNanobeamFEM.m:420-427 cross-checks pec_inds against
% COMSOL's DEFAULT PEC selection and hard-errors on any index that is not in
% it. Every face in these selections is exterior, so that check passes; the
% check is one-directional (extra defaults are fine), which is what lets the
% shield's own outer faces keep the default PEC.
if useOpt
    % x = 0: beam cross-section plus the air's half-disc face.
    cylXsymSel = beamgeom.create('cylXsymSel', 'BoxSelection');
    cylXsymSel.set('xmin', -d2).set('xmax', d2);
    cylXsymSel.set('ymin', -d2).set('ymax', airrad + d2);
    cylXsymSel.set('zmin', zLoCyl - d2).set('zmax', airrad + d2);
    cylXsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    P.bndSel.cylXsym = readSel(model, P.geomname, 'cylXsymSel');

    % y = 0: the whole plane out to xS1, so the beam's, pad's and shield's
    % y = 0 faces get the optical symmetry BC too and not just the air's. The
    % holes' own y = 0 faces are in here as well (each hole straddles y = 0 in
    % the full device, so in the half model its y = 0 cut is an exterior face
    % on the mirror plane). xS1 rather than xP1 keeps the PML out: the PML
    % domains are not in the emw selection.
    cylYsymSel = beamgeom.create('cylYsymSel', 'BoxSelection');
    cylYsymSel.set('xmin', -d2).set('xmax', xS1 + d2);
    cylYsymSel.set('ymin', -d2).set('ymax', d2);
    cylYsymSel.set('zmin', zLoCyl - d2).set('zmax', airrad + d2);
    cylYsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    P.bndSel.cylYsym = readSel(model, P.geomname, 'cylYsymSel');

    % z = 0: exactly structZsymSel, i.e. solid AND air faces. Built above
    % inside the symZOn branch, so reuse it rather than emitting a duplicate
    % box. Only defined when symZOn, which readCrossShieldParams guarantees
    % whenever abs(P.oevenz) = 1.
    if symZOn
        P.bndSel.cylZsym = readSel(model, P.geomname, 'structZsymSel');
    end

    % x = airCylLen end cap, EXTERIOR faces only.
    %
    % cylXendAllSel catches everything in the plane: the cap itself (vacuum
    % outside - a genuine scattering face) and the imprint of the cap on the
    % beam (diamond on BOTH sides - an interior face, and a Scattering BC on an
    % interior face is an absorber buried inside the cavity mirror).
    % cylXendBeamSel catches only the imprint, because the beam cross-section
    % lies inside y <= wid/2, |z| <= thi/2 while every face of the cap has
    % vertices out on the circle at radius airrad.
    %
    % This is also where BuildNanobeamFEM.m:443-444's bug lives - cylXendSel
    % there sets ymin/ymax twice and never sets zmin/zmax, so the z bounds fall
    % back to the defaults. Not replicated: both boxes below set all six.
    cylXendAllSel = beamgeom.create('cylXendAllSel', 'BoxSelection');
    cylXendAllSel.set('xmin', xC1 - d2).set('xmax', xC1 + d2);
    cylXendAllSel.set('ymin', -d2).set('ymax', airrad + d2);
    cylXendAllSel.set('zmin', zLoCyl - d2).set('zmax', airrad + d2);
    cylXendAllSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;

    cylXendBeamSel = beamgeom.create('cylXendBeamSel', 'BoxSelection');
    cylXendBeamSel.set('xmin', xC1 - d2).set('xmax', xC1 + d2);
    cylXendBeamSel.set('ymin', -d2).set('ymax', wid/2 + d2);
    cylXendBeamSel.set('zmin', zLoEff - d2).set('zmax', thi/2 + d2);
    cylXendBeamSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;

    cylXendSel = beamgeom.create('cylXendSel', 'DifferenceSelection');
    cylXendSel.set('entitydim', 2);
    cylXendSel.set('add', 'cylXendAllSel').set('subtract', 'cylXendBeamSel');
    beamgeom.runCurrent;
    P.bndSel.cylXend = readSel(model, P.geomname, 'cylXendSel');

    nXendAll  = numel(readSel(model, P.geomname, 'cylXendAllSel'));
    nXendBeam = numel(readSel(model, P.geomname, 'cylXendBeamSel'));
    if nXendBeam == 0
        warning('BuildNanobeamCrossShieldFEM:noEndCapImprint', ...
            ['The x = %.1f nm plane contains no face inside the beam ' ...
             'cross-section, so the air end cap apparently did not imprint on ' ...
             'the beam. Either Compose dropped the interior boundary - in ' ...
             'which case beam and air have MERGED into one domain and the ' ...
             'materials are wrong - or the cap landed in a hole. Check the ' ...
             'geometry in the GUI.'], xC1*1e9);
    end
    disp(['Air end cap at x = ', num2str(xC1*1e9,'%.0f'), ' nm: ', ...
          num2str(nXendAll), ' faces in plane, ', num2str(nXendBeam), ...
          ' interior (removed), ', num2str(numel(P.bndSel.cylXend)), ...
          ' scattering']);

    % Curved outer surface, found with the 45 degree 'intersects' probe of
    % BuildNanobeamFEM.m:470-477. The point (y,z) = airrad/sqrt(2)*(1,1) is on
    % the cylinder whether or not the z < 0 half has been cut away, so the same
    % probe works for both symZOn cases. The non-rect second probe
    % (cylCurv2Sel) is not ported - v1 is rect only.
    cylCurvSel = beamgeom.create('cylCurvSel', 'BoxSelection');
    cylCurvSel.set('xmin', d2).set('xmax', 2*d2);
    cylCurvSel.set('ymin', airrad/sqrt(2) - 10*d2).set('ymax', airrad/sqrt(2) + 10*d2);
    cylCurvSel.set('zmin', airrad/sqrt(2) - 10*d2).set('zmax', airrad/sqrt(2) + 10*d2);
    cylCurvSel.set('entitydim', 2).set('condition', 'intersects');
    beamgeom.runCurrent;
    P.bndSel.cylCurv = readSel(model, P.geomname, 'cylCurvSel');
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
if useOpt
    % Every one of these is read unconditionally or near-unconditionally by
    % SetupNanobeamFEM's optical block; cylXend and cylCurv are the two that
    % are ALWAYS read (they make up bnds.scat_inds at :390), so an empty one
    % there means the cavity has no radiating boundary at all and the optical Q
    % comes out meaninglessly high.
    assertNonEmptySel('P.domSel.cyl',     P.domSel.cyl);
    assertNonEmptySel('P.bndSel.cylXsym', P.bndSel.cylXsym);
    assertNonEmptySel('P.bndSel.cylYsym', P.bndSel.cylYsym);
    assertNonEmptySel('P.bndSel.cylXend', P.bndSel.cylXend);
    assertNonEmptySel('P.bndSel.cylCurv', P.bndSel.cylCurv);
    if symZOn
        assertNonEmptySel('P.bndSel.cylZsym', P.bndSel.cylZsym);
    end

    % The mechanical z = 0 selection is the optical one minus the air faces, so
    % it must be a strict subset and must not have lost everything.
    if symZOn
        if ~isempty(setdiff(P.bndSel.beamZsym, P.bndSel.cylZsym))
            warning('BuildNanobeamCrossShieldFEM:ZsymNotSubset', ...
                ['P.bndSel.beamZsym contains boundaries [%s] that are not in ' ...
                 'P.bndSel.cylZsym. The two are built by differencing the ' ...
                 'same box, so this cannot happen unless a selection feature ' ...
                 'was renamed.'], ...
                num2str(setdiff(P.bndSel.beamZsym, P.bndSel.cylZsym)));
        end
        nAirZ0 = numel(P.bndSel.cylZsym) - numel(P.bndSel.beamZsym);
        if nAirZ0 <= 0
            warning('BuildNanobeamCrossShieldFEM:noAirZsymRemoved', ...
                ['The z = 0 air-face difference removed %d boundaries. The ' ...
                 'air cylinder must have z = 0 faces (the outer quarter-disc ' ...
                 'region plus one per enclosed hole), so 0 means the y ' ...
                 'discrimination in airZ0OutSel/airZ0HoleSel missed them and ' ...
                 'the Solid Mechanics symmetry BC is about to be hung on air.'], ...
                nAirZ0);
        else
            disp(['z = 0 plane: ', num2str(numel(P.bndSel.cylZsym)), ...
                  ' faces total, ', num2str(nAirZ0), ' air (optical only), ', ...
                  num2str(numel(P.bndSel.beamZsym)), ' solid (mechanical)']);
        end
    end

    % The air must not have been swept into the solid domain selections. This
    % is the cheap version of the z-extent argument at the top of the domain
    % selection block, checked against the geometry COMSOL actually built.
    if ~isempty(intersect(P.domSel.cyl, P.domSel.beam))
        warning('BuildNanobeamCrossShieldFEM:airInStructSel', ...
            ['Domain(s) [%s] are in BOTH P.domSel.cyl and P.domSel.beam. The ' ...
             'air would then be given the beam material and Solid Mechanics ' ...
             'would be solved in it.'], ...
            num2str(intersect(P.domSel.cyl, P.domSel.beam)));
    end
    if usePML && ~isempty(intersect(P.domSel.cyl, P.domSel.PML))
        warning('BuildNanobeamCrossShieldFEM:airInPMLSel', ...
            'Domain(s) [%s] are in both P.domSel.cyl and P.domSel.PML.', ...
            num2str(intersect(P.domSel.cyl, P.domSel.PML)));
    end
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

% --- optical fields ------------------------------------------------------
% Nothing below this line runs unless P.solveOpt = 1, so the P.solveOpt = 0
% behaviour of this function is unchanged.
if isfield(P,'solveOpt') && P.solveOpt
    oRequired = {'airrad','lambda','nbeam','oevenx','oeveny','oevenz'};
    oMissing  = oRequired(~isfield(P, oRequired));
    if ~isempty(oMissing)
        error('BuildNanobeamCrossShieldFEM:missingOpticalFields', ...
            ['P.solveOpt = 1 needs P field(s): %s. The house convention for ' ...
             'the radius is P.airrad = 2*P.lambda + P.w/2.'], ...
            strjoin(oMissing, ', '));
    end
    validateattributes(P.airrad, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.airrad');
    validateattributes(P.lambda, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.lambda');
    validateattributes(P.nbeam,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.nbeam');

    xB1 = P.beamLenHalf;
    xS0 = xB1 + P.shieldPadLen;

    % Default cylinder extent: as far into the hole mirror as it can reach
    % while still leaving one free-space wavelength of vacuum before the shield
    % front face, and never past the beam end (a cap buried inside the clamp
    % pad would be cut to the shape of the pad cross-section instead of the
    % beam's, and cylXendBeamSel is sized for the beam).
    if ~isfield(P,'airCylLen') || isempty(P.airCylLen)
        P.airCylLen = min(xB1, xS0 - P.lambda);
    end
    validateattributes(P.airCylLen, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.airCylLen');

    % THE guard that replaces the old blanket P.solveOpt rejection. A hard
    % error, not a warning: an overlapping cylinder silently turns the phononic
    % shield into a lump of diamond in an air bath, and the Q it then reports
    % looks entirely plausible.
    if P.airCylLen >= xS0
        error('BuildNanobeamCrossShieldFEM:airCylReachesShield', ...
            ['P.airCylLen = %.1f nm reaches the shield front face at ' ...
             'xS0 = %.1f nm (= P.beamLenHalf %.1f + P.shieldPadLen %.1f). ' ...
             'The air cylinder would then enclose part of the shield, which ' ...
             'is exactly the failure the old P.solveOpt guard existed to ' ...
             'prevent. Reduce P.airCylLen below %.1f nm.'], ...
            P.airCylLen*1e9, xS0*1e9, xB1*1e9, P.shieldPadLen*1e9, xS0*1e9);
    end

    % z = 0 exists as a geometric cut, or it does not; there is one cut shared
    % by both physics. SetupNanobeamFEM hangs a mechanical symmetry BC on
    % P.bndSel.beamZsym whenever abs(P.mevenz) = 1 and an optical PMC/PEC on
    % P.bndSel.cylZsym whenever abs(P.oevenz) = 1, and both selections only
    % exist when the cut is made.
    if (abs(P.oevenz) > 0) ~= (abs(P.mevenz) > 0)
        error('BuildNanobeamCrossShieldFEM:evenzMismatch', ...
            ['P.oevenz = %g and P.mevenz = %g disagree about whether z = 0 is ' ...
             'a mirror plane. There is only one z = 0 cut in the geometry, so ' ...
             'either both must be +/-1 (eighth model) or both must be 0 ' ...
             '(quarter model). As written, one physics would get a symmetry ' ...
             'BC on a selection that resolves to nothing.'], ...
            P.oevenz, P.mevenz);
    end

    % The air's z = 0 faces are separated from the solid ones by y alone (see
    % the airZ0OutSel/airZ0HoleSel comment). That only works while the holes
    % stay clear of the beam side wall, which is also the condition for the
    % holes not to break through and make the beam a comb.
    maxHoleHalfWid = max(P.geomHalf(:,2))/2;
    if maxHoleHalfWid >= P.w/2 - 10e-9
        error('BuildNanobeamCrossShieldFEM:holesReachSideWall', ...
            ['max(hy)/2 = %.1f nm is within 10 nm of w/2 = %.1f nm. The ' ...
             'z = 0 air faces inside the holes can then no longer be told ' ...
             'apart from the beam''s own z = 0 face by a y bound, so ' ...
             'P.bndSel.beamZsym would pick up air boundaries. Narrow P.hy or ' ...
             'widen P.w.'], maxHoleHalfWid*1e9, P.w/2*1e9);
    end
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

function assertCylClears(P, xC1, airrad, zLoCyl, xB1, xS0, xS1, yS1, ...
                         xP1, yP1, zLoEff, thi, usePML)
%ASSERTCYLCLEARS Hard error if the air cylinder envelope touches shield or PML.
%
% This is the replacement for the blanket "P.solveOpt is not supported" guard,
% and it is a hard error for the same reason that guard was: a cylinder that
% overlaps the shield does not fail, it just quietly dissolves the phononic
% bandgap into an air bath and reports a plausible-looking Q.
%
% The test is on axis-aligned bounding boxes, which is conservative for a
% cylinder (the AABB is the circumscribing box) and therefore safe. Each
% candidate is reported with the arithmetic, because for this geometry the y
% and z extents genuinely DO overlap and x is the only thing keeping the
% cylinder out of the shield - a reader needs to see which term saved them.

cylBox = [0, xC1; 0, airrad; zLoCyl, airrad];

boxes = {'cross shield', [xS0, xS1; 0, yS1; zLoEff, thi/2]};
if usePML
    boxes = [boxes, { ...
        'PML +x arm',  [xS1, xP1; 0,   yS1; zLoEff, thi/2], ...
        'PML +y arm',  [xB1, xS1; yS1, yP1; zLoEff, thi/2], ...
        'PML corner',  [xS1, xP1; yS1, yP1; zLoEff, thi/2]}];
end

axisName = {'x','y','z'};
for k = 1:2:numel(boxes)
    name  = boxes{k};
    other = boxes{k+1};

    % Positive separation on ANY axis means the boxes are disjoint. Collect the
    % gaps so the message can name the axis that is doing the work.
    gap = zeros(1,3);
    for ax = 1:3
        gap(ax) = max(other(ax,1) - cylBox(ax,2), cylBox(ax,1) - other(ax,2));
    end

    if all(gap <= 0)
        [~, best] = max(gap);
        error('BuildNanobeamCrossShieldFEM:airCylOverlapsStructure', ...
            ['The air cylinder envelope overlaps the %s.\n' ...
             '  cylinder  x [%.1f %.1f]  y [%.1f %.1f]  z [%.1f %.1f] nm\n' ...
             '  %-12s x [%.1f %.1f]  y [%.1f %.1f]  z [%.1f %.1f] nm\n' ...
             '  closest separation is on %s and is %.1f nm (needs > 0).\n' ...
             'Reduce P.airCylLen (now %.1f nm) or P.airrad (now %.1f nm). ' ...
             'This is the failure the old blanket P.solveOpt guard prevented.'], ...
            name, ...
            cylBox(1,1)*1e9, cylBox(1,2)*1e9, cylBox(2,1)*1e9, cylBox(2,2)*1e9, ...
            cylBox(3,1)*1e9, cylBox(3,2)*1e9, ...
            name, ...
            other(1,1)*1e9, other(1,2)*1e9, other(2,1)*1e9, other(2,2)*1e9, ...
            other(3,1)*1e9, other(3,2)*1e9, ...
            axisName{best}, gap(best)*1e9, xC1*1e9, airrad*1e9);
    end

    sepStr = '';
    for ax = 1:3
        if gap(ax) > 0
            sepStr = [sepStr, sprintf(' %s by %.1f nm', axisName{ax}, gap(ax)*1e9)]; %#ok<AGROW>
        end
    end
    disp(['Air cylinder clears the ', name, ':', sepStr]);
end

% Informational, not fatal: the cylinder is allowed to be taller in y than the
% shield is, because the x separation is what matters, but it is worth knowing.
% It is normal and expected for a deliberately small shield (P.nShieldY = 1 or
% 2, as in the optics-only bring-up stage); it is a warning sign for a
% production shield, because it means the air is wider than the device.
if airrad > yS1
    warning('BuildNanobeamCrossShieldFEM:airCylTallerThanShield', ...
        ['P.airrad = %.1f nm exceeds the shield height yS1 = ' ...
         'P.nShieldY*P.aShield = %.1f nm. The geometry is still valid - the ' ...
         'cylinder and the shield are separated in x - but the air now ' ...
         'extends in y past the whole shield/PML frame, so nothing is meshed ' ...
         'alongside the frame''s inboard flank for x < P.airCylLen. Expected ' ...
         'if you have shrunk P.nShieldY for a bring-up run.'], ...
        airrad*1e9, yS1*1e9);
end

% Cutting a hole is not fatal but it breaks the cylXend difference: the air
% inside a cut hole has a cross-section that lies entirely within the beam
% box, so cylXendBeamSel would subtract it and the cap would lose a piece.
xh0 = P.geomHalf(:,3) - P.geomHalf(:,1)/2;
xh1 = P.geomHalf(:,3) + P.geomHalf(:,1)/2;
cut = find(xh0 < xC1 & xh1 > xC1);
if ~isempty(cut)
    warning('BuildNanobeamCrossShieldFEM:airCylCutsHole', ...
        ['P.airCylLen = %.1f nm falls inside hole %s (x in [%.1f, %.1f] nm). ' ...
         'The end cap then includes that hole''s cross-section, which lies ' ...
         'inside the beam bounding box and will be differenced OUT of ' ...
         'P.bndSel.cylXend - leaving part of the cap with no scattering BC. ' ...
         'Move P.airCylLen into the solid between two holes.'], ...
        xC1*1e9, mat2str(cut'), xh0(cut(1))*1e9, xh1(cut(1))*1e9);
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

function auditCrossShield(P, wid, thi, usePML, useOpt)
%AUDITCROSSSHIELD Pure-MATLAB, warn-only checks run before anything is built.
%
% Nothing here changes the geometry. The point is that every compromise in this
% design is numeric and knowable in advance, so none of them should be silent.
%
% useOpt gates the optical block at the end; the checks above it are unchanged.

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

% --- optical: cylinder placement and the combined mesh cost --------------
if useOpt
    ac  = P.airCyl;
    lam = P.lambda;

    disp(['Air cylinder: r = ', num2str(ac.rad*1e9,'%.0f'), ' nm, x in [0, ', ...
          num2str(ac.len*1e9,'%.0f'), '] nm, ', num2str(ac.gapToxS0*1e9,'%.0f'), ...
          ' nm (', num2str(ac.gapToxS0/lam,'%.2f'), ' free-space wavelengths) ', ...
          'of bare vacuum before the shield, ', num2str(ac.nHolesIn), ' of ', ...
          num2str(size(P.geomHalf,1)), ' half-beam holes enclosed']);

    % The mode has to be dead by the time the air stops, because beyond that
    % the beam is clad in the EM interface's default PEC - a perfect mirror
    % where the real device has vacuum. ndef holes are the defect, so
    % nHolesIn - ndef is the number of mirror periods actually inside the air.
    nMirrorIn = ac.nHolesIn - localGetField(P,'ndef',0);
    if nMirrorIn < 6
        warning('BuildNanobeamCrossShieldFEM:airCylTooShort', ...
            ['Only %d mirror period(s) lie inside the air cylinder ' ...
             '(%d holes enclosed, %g of them defect holes). The beam beyond ' ...
             'x = %.1f nm is clad in the default PEC of the EM interface, so ' ...
             'with this few mirror periods the cavity is terminated by a ' ...
             'perfect mirror and the optical Q is a fiction. Increase ' ...
             'P.airCylLen, or add holes.'], ...
            nMirrorIn, ac.nHolesIn, localGetField(P,'ndef',0), ac.len*1e9);
    end

    % Combined mesh cost. Same crude tet volume h^3/6 as nElemEst above, which
    % understates the count by about 1.4x against a real Delaunay fill - good
    % enough to tell "fits" from "does not fit", not good enough to tune with.
    %
    % The two physics do NOT share degrees of freedom: SetupNanobeamFEM puts
    % emw on P.domSel.beam + P.domSel.cyl and smech on P.domSel.beam +
    % P.domSel.PML, so the air never contributes a mechanical DOF and the PML
    % never contributes an optical one. They only share the mesh.
    if abs(P.mevenz) > 0
        zFrac = 0.5;        % eighth model: air is a quarter cylinder
    else
        zFrac = 1.0;        % quarter model: air is a half cylinder
    end
    hOpt   = lam/5;                                  % SolveNanobeamFEM.m:131
    vAir   = zFrac*pi*ac.rad^2/2 * ac.len;
    nTetO  = vAir/(hOpt^3/6);
    dofO   = 20*nTetO;                               % 2nd order curl elements

    if usePML
        vPML = (P.PMLLen*P.shield.yS1 + ...
                (P.shield.xS1 - P.shield.xB1)*P.PMLLenY + ...
                P.PMLLen*P.PMLLenY) * thi*zFrac;
        if isfield(P,'PMLmeshDiv') && ~isempty(P.PMLmeshDiv) && P.PMLmeshDiv > 0
            hPML = P.PMLLen/P.PMLmeshDiv;
        else
            hPML = P.shieldHmax;
        end
        nTetPML = vPML/(hPML^3/6);
    else
        nTetPML = 0;
    end
    nTetM = nElemEst + nTetPML;
    dofM  = 4.2*nTetM;          % 3 components x ~1.4 quadratic nodes per tet

    disp(['Combined mesh estimate: optical ', num2str(nTetO,'%.2e'), ' tets / ', ...
          num2str(dofO,'%.2e'), ' DOF at hmax = ', num2str(hOpt*1e9,'%.0f'), ...
          ' nm; mechanical ', num2str(nTetM,'%.2e'), ' tets / ', ...
          num2str(dofM,'%.2e'), ' DOF at shield hmax = ', ...
          num2str(P.shieldHmax*1e9,'%.1f'), ' nm. P.max_dof = ', ...
          num2str(localGetField(P,'max_dof',NaN),'%.1e')]);

    maxDof = localGetField(P,'max_dof',Inf);
    if dofM > maxDof
        warning('BuildNanobeamCrossShieldFEM:mechDofOverBudget', ...
            ['The mechanical estimate (%.2e DOF) is already over P.max_dof ' ...
             '(%.1e) BEFORE optics - the shield hmax is doing it, not the air ' ...
             '(the air contributes 0 mechanical DOF). SolveNanobeamFEM''s ' ...
             'mAdjMesh loop cannot fix this: it raises hauto, but ' ...
             'applyMeshOverrides re-applies hmax = P.shieldHmax on the shield ' ...
             'every pass, so the dominant term never moves. Levers, in order ' ...
             'of preference: reduce P.nShieldX/P.nShieldY for the combined ' ...
             'run, relax P.shieldHmax towards - but not all the way to - ' ...
             '(a-h)/2 = %.1f nm, or raise P.max_dof.'], ...
            dofM, maxDof, wall/2*1e9);
    end
    if dofO > maxDof
        warning('BuildNanobeamCrossShieldFEM:optDofOverBudget', ...
            ['The optical estimate (%.2e DOF) is over P.max_dof (%.1e). ' ...
             'P.oAdjMesh will coarsen the WHOLE mesh, air included, until it ' ...
             'fits. Reduce P.airrad or P.airCylLen.'], dofO, maxDof);
    end
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
