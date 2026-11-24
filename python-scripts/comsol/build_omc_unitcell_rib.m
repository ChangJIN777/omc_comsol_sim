function model = build_omc_unitcell_rib(outFile, xsect)
% BUILD_OMC_UNITCELL_RIB  Build the parametric rib unit-cell template
% (omc_unitcell_rib.mph) via COMSOL LiveLink for MATLAB.
%
%   model = build_omc_unitcell_rib()         % saves comsol/omc_unitcell_rib.mph
%   model = build_omc_unitcell_rib(outFile)  % saves to outFile
%   model = build_omc_unitcell_rib('')       % build in memory only, no save
%   model = build_omc_unitcell_rib('wedge')  % saves comsol/omc_unitcell_rib_wedge.mph
%   model = build_omc_unitcell_rib(outFile, 'wedge')   % wedge, saved to outFile
%   model = build_omc_unitcell_rib('', 'wedge')        % wedge, in memory only
%
% XSECT (cross-section; python-scripts config key `xsect`):
%   'rect'   (default) the slab only: vertical walls, flat underside.
%   'wedge'  the slab plus the underside left by an angled etch: the "roof"
%            of the FILLETED outline. Under a point p of the cell the wedge
%            is d(p)*tan(wedge_angle) deep, where d(p) is the in-plane
%            distance from p to the filleted outline (neighbouring cells
%            included) and wedge_angle = 30 deg. On a straight strip of
%            width b this is the wedge of the hole template
%            omc_unitcell_iso.mph (its Polygon 1: (-t/2, +-w/2),
%            (-t/2 - w/2/sqrt(3), 0)): an isosceles triangle on the bottom
%            face of the slab, base angles 30 deg, apex b/(2*sqrt(3)) below.
%            With u = t/2 - |z|, v = s/2 - |y|, f = w/2 - |y|:
%              d = min( hypot(max(u+r1,0), max(v+r1,0)) - r1 , f )
%              d = r2 - hypot(r2-u, r2-f)   in the tip corners (0<=u<r2, f<r2)
%            (python-scripts/src/geometry.py: rib_roof_distance). It is built
%            EXACTLY from planes and cones, in zones that meet on planes
%            only, where their profiles agree (no tangent intersections):
%              S   |z| >= t/2+r1            spine prism, b = s
%              J   |z| <= t/2+r1, y <= s/2+r1   box minus a cone frustum about
%                  the root-arc centre (radius r1 at the slab, widening by
%                  1/tan(wedge_angle) per depth): d = rho - r1. Covers the
%                  root fillet AND the whole junction, the deepest place.
%              A   y >= s/2+r1              rib prism along y, b = t
%              all three are then cut by the tip plane d = w/2 - y (the hip),
%              the two tip-corner boxes are removed, and
%              C   quarter cones (radius r2 at the slab, apex r2*tan(angle)
%                  below the tip-arc centre) are put in their place.
%            Sequence: outline -> fillets (wp1) -> slab; wedge from the same
%            r1, r2 arcs -> union with the slab. The wedge has no
%            mid-thickness mirror plane. NOT BUILT OR SOLVED YET: see the
%            checks listed at the wedge block below.
%   A first argument 'rect' or 'wedge' is taken as XSECT, not as a file name.
%
% Rib counterpart of omc_unitcell.mph. The node structure (tags, labels,
% physics, boundary conditions, material, mesh approach, studies) is copied
% from omc_unitcell.mph; the differences are the geometry, 20 (not 12)
% eigenfrequencies in the two mechanical studies, and ANISOTROPIC (cubic)
% diamond with C11/C12/C44 as global parameters. The crystal is rotated by
% the global parameter rxt about the thickness axis x: beam along <100> for
% rxt = 0 (default), along <110> for rxt = 45[deg].
%
% AIR: a half-disk (half-cylinder, y >= 0) around the beam with its axis
% along z, as Work Plane 2 / Revolve 1 of omc_unitcell_iso.mph (radius 2 um
% about the mid-depth of its wedge) and the air cylinder of
% omc-comsol-chang/buildRibUnitCell_LN.m (radius 5a). omc_unitcell.mph has a
% box instead (5*t by 1.5*w). The axis is at the mid-depth of the beam
% (x = air_xc), the radius air_r reaches air_gap beyond the farthest point
% of the beam, and the curved surface carries the scattering condition
% sctr1. Templates built before this change have the box and no parameter
% air_r: its faces x = +-2.5*th were about 2 decay lengths from the slab,
% inside the evanescent tail of the guided modes.
% air_gap is a global parameter and every air size derives from it (through
% air_r), so the air can be resized without a rebuild: in the GUI (Global
% Definitions > Parameters 1), with model.param.set('air_gap', ...), or in
% a batch solve with -p air_gap=2.5[um] (python-scripts/scripts/
% run_comsol_batch.sh).
%
% GEOMETRY: the 'rib' cell of omc-comsol-chang/buildRibUnitCell_LN.m with
% vertical walls, i.e. the per-cell recipe of addRibCell.m: a continuous
% spine of width s carrying one rib per cell, centred on the cell, t wide
% along the beam and reaching |y| = w/2 = L from the beam axis; r2 rounds
% the two convex rib-tip corners and r1 the two concave rib/spine corners
% (per half). Two changes from addRibCell.m:
%   * every dimension is a parameter EXPRESSION, so the saved .mph is tunable
%   * the frame is the template's: beam axis z, width y, thickness x
%     (buildRibUnitCell_LN uses x, y, z), so Floquet k = (0,0,kF) is unchanged
% Names follow the 'rib' P struct (omc-comsol-chang/rib_README.md).
%
% Global parameters (all referenced by name in the geometry):
%   a   period along z (cell width)         th  slab thickness (x)
%   w   rib tip-to-tip width (y), = 2L;     t   rib width along z (beam axis)
%       the model draws the y >= 0 half     s   spine width (y), full width
%   r1  concave root fillet (rib/spine)     r2  convex rib-tip fillet
%   kF  Floquet wavevector along z (swept by the driver)
%   n_diamond    diamond refractive index, 2.386 (telecom; the value of
%                python-scripts/configs/materials.yaml and of P.nbeam in
%                omc-comsol-chang). relpermittivity = n_diamond^2.
%   C11, C12, C44  cubic stiffness of diamond (configs/materials.yaml)
%   rxt          crystal rotation about the thickness axis x: the beam is at
%                rxt from <100> (P.rxtal of omc-comsol-chang). ONLY multiples
%                of 45 deg: any other angle couples normal and shear strain
%                across the y = 0 plane and breaks the y mirror of the half
%                model (sym1 / as1).
%   Cani         C11-C12-2*C44, the cubic anisotropy constant
%   dCrot        Cani/2*sin(2*rxt)^2, in-plane stiffness change of the
%                rotated crystal
%   gCrot        Cani/4*sin(4*rxt), normal/shear coupling of the rotated
%                crystal; 0 at multiples of 45 deg
%   wedge_angle  base angle of the wedge, 30[deg] ('wedge' only)
%   air_gap      air between the beam and the scattering boundary, 1.55[um]
%   air_xc       x of the air half-disk axis: 0 ('rect'), -wd/2 ('wedge')
%   air_r        air half-disk radius, sqrt((w/2)^2+(th/2-air_xc)^2)+air_gap
% r1 and r2 must be > 0 (COMSOL's Fillet rejects radius 0); use ~1 nm for a
% near-sharp corner. The geometry needs t < a, s < w, 2*r2 < t,
% t + 2*r1 < a and r1 + r2 < (w - s)/2 (a straight rib flank).
%
% Half cell seen from +x (y up, z right):
%
%   y = w/2              +---t---+          r2: both rib-tip corners
%                        |       |
%                        |  rib  |          r1: both rib/spine corners
%                        |       |
%   y = s/2      +-------+       +-------+
%                |        spine          |
%   y = 0        +-----------------------+   <- y symmetry plane
%              -a/2    -t/2     t/2     a/2     (z)
%
% Studies (as omc_unitcell.mph). The model is the y >= 0 half; the y = 0
% plane carries a symmetry or an antisymmetry condition, one per study:
%   std1 'mech sym'      Solid, sym1 (y-even modes)      -> dataset dset1
%   std3 'mech antisym'  Solid, as1  (y-odd modes)       -> dataset dset2
%   std2 'opt TE'        EWFD, pec2 on y = 0             -> dataset dset3
%   std4 'opt TM'        EWFD, pmc1 on y = 0             -> dataset dset4
% Solver sequences are created here, in that order, so the dataset of each
% study is fixed (omc_unitcell.mph has none and gets dset1 from whichever
% study runs first).
%
% Deliberate differences from buildRibUnitCell_LN.m:
%   * vertical walls only (its sidewall angle d, addRibPrisms.m, is not
%     carried over).
%   * no thickness-symmetry cut (zSymCondition). omc_unitcell.mph has none.
%   * no fixed-length L(a) rule: w is a plain parameter (the builder's
%     explicit-P.w path, as test_rib_unitCell_explicitW.m).
%   * selections are parametric boxes that follow the drawn cell
%     (z = -a/2..a/2, y = 0..w/2), and the geometry is finalized with run.
%
% Requires: COMSOL 6.2+ with LiveLink for MATLAB connected (mphstart), the
% Structural Mechanics (or MEMS / Acoustics) and Wave Optics modules.
%
% See also OMC_MECHANICAL_LIVELINK_RIB, OMC_MECHANICAL_LIVELINK.

import com.comsol.model.*
import com.comsol.model.util.*

% Arguments: (outFile), (xsect), (outFile, xsect).
if nargin >= 1 && (ischar(outFile) || isstring(outFile)) ...
        && any(strcmp(char(outFile), {'rect', 'wedge'})) && nargin < 2
    xsect = char(outFile);
    haveFile = false;
else
    haveFile = nargin >= 1;
    if nargin < 2
        xsect = 'rect';
    end
end
xsect = char(xsect);
if ~any(strcmp(xsect, {'rect', 'wedge'}))
    error('build_omc_unitcell_rib:xsect', ...
        'xsect must be ''rect'' or ''wedge'', got ''%s''.', xsect);
end
isWedge = strcmp(xsect, 'wedge');
if isWedge
    defaultName = 'omc_unitcell_rib_wedge.mph';
else
    defaultName = 'omc_unitcell_rib.mph';
end
if ~haveFile
    outFile = fullfile(fileparts(mfilename('fullpath')), defaultName);
end
if ~isempty(outFile) && exist(outFile, 'file')
    error('build_omc_unitcell_rib:fileExists', ...
        '%s already exists. Delete it or pass another file name.', outFile);
end

nEigMech = 20;      % eigenfrequencies per mechanical study (template: 12)

%% Model and global parameters
model = ModelUtil.create(char(ModelUtil.uniquetag('Model')));
model.label(defaultName);

% Defaults: omc-comsol-chang/test_rib_unitCell_explicitW.m (diamond, d = 0).
params = { ...
    'a',  '570[nm]',  'lattice constant = cell width (period along z)'; ...
    'w',  '1747[nm]', 'rib tip-to-tip width (y), = 2L; the y >= 0 half is drawn'; ...
    'th', '283[nm]',  'slab thickness (x)'; ...
    't',  '391[nm]',  'rib width along the beam (z)'; ...
    's',  '254[nm]',  'spine width (y), full width'; ...
    'r1', '20[nm]',   'concave rib/spine (root) fillet radius, > 0'; ...
    'r2', '20[nm]',   'convex rib-tip fillet radius, > 0'; ...
    'kF', 'pi/a',     'Floquet wavevector along z (swept by driver)'; ...
    'n_diamond', '2.386', 'diamond refractive index, telecom (configs/materials.yaml)'; ...
    'C11', '1076[GPa]', 'diamond cubic stiffness C11 (configs/materials.yaml)'; ...
    'C12', '125[GPa]',  'diamond cubic stiffness C12 (configs/materials.yaml)'; ...
    'C44', '577[GPa]',  'diamond cubic stiffness C44 (configs/materials.yaml)'; ...
    'rxt', '0[deg]', ...
        'crystal rotation about thickness x: beam at rxt from <100>; 0 or 45[deg] only (half model)'; ...
    'Cani', 'C11-C12-2*C44', 'cubic anisotropy constant (0 for an isotropic solid)'; ...
    'dCrot', 'Cani/2*sin(2*rxt)^2', 'in-plane stiffness change of the rotated crystal'; ...
    'gCrot', 'Cani/4*sin(4*rxt)', ...
        'normal/shear coupling of the rotated crystal; 0 at multiples of 45[deg]'};
for i = 1:size(params, 1)
    model.param.set(params{i, 1}, params{i, 2}, params{i, 3});
end
if isWedge
    model.param.set('wedge_angle', '30[deg]', ...
        'wedge base angle: apex b/2*tan(wedge_angle) below a strip of width b');
    model.param.set('wdj', 'tan(wedge_angle)*(sqrt((t/2+r1)^2+(s/2+r1)^2)-r1)', ...
        'wedge depth under the junction centre (root cone, before the tip cut)');
    model.param.set('wd', 'min(wdj,tan(wedge_angle)*w/2)', ...
        'deepest point of the wedge below the slab');
end
% Air half-disk (built after the beam, below). air_xc is the mid-depth of
% the beam, so the air above the top face and below the wedge apex is equally
% thick. air_r reaches air_gap beyond the farthest point of the beam from
% that axis, the top corner of the rib tip (x = th/2, y = w/2).
if isWedge
    airXc = '-wd/2';
else
    airXc = '0[nm]';
end
model.param.set('air_gap', '1.55[um]', ...
    'air between the beam and the scattering boundary (one 1550 nm wavelength)');
model.param.set('air_xc', airXc, ...
    'x of the air half-disk axis: mid-depth of the beam');
model.param.set('air_r', 'sqrt((w/2)^2+(th/2-air_xc)^2)+air_gap', ...
    'radius of the air half-disk (axis along z, y >= 0)');

%% Geometry (buildRibUnitCell_LN / addRibCell recipe, template frame)
comp = model.component.create('comp1', true);
geom = comp.geom.create('geom1', 3);
% COMSOL kernel, set explicitly: without this the geometry takes the user's
% preference (Geometry > Geometry representation in new models). With the
% CAD kernel ('cadps', Parasolid) every build AND every later solve of the
% saved .mph needs a CAD Import Module seat (CADIMPORT / CADIMPORTBATCH).
% omc_unitcell.mph and omc_unitcell_iso.mph are 'comsol' as well.
geom.geomRep('comsol');
geom.lengthUnit('m');

% Work plane x = -th/2. Its local axes are (y, z); the normal is +x.
wp = geom.create('wp1', 'WorkPlane');
wp.set('unite', true);
wp.set('planetype', 'quick');
wp.set('quickplane', 'yz');
wp.set('quickx', '-th/2');

% addRibCell.m feature chain, y >= 0 half, centred on z = 0:
%   r1_u  rib [w/2 x t]   (y = 0..w/2,  z = -t/2..t/2)
%   r3_u  half-spine      (y = 0..s/2,  z = -a/2..a/2)
rib = wp.geom.create('r1_u', 'Rectangle');
rib.set('base', 'corner');
rib.set('pos',  {'0', '-t/2'});
rib.set('size', {'w/2', 't'});
spine = wp.geom.create('r3_u', 'Rectangle');
spine.set('base', 'corner');
spine.set('pos',  {'0', '-a/2'});
spine.set('size', {'s/2', 'a'});

uni = wp.geom.create('union_u', 'Union');
uni.set('intbnd', false);
uni.selection('input').set({'r1_u', 'r3_u'});

% Fillets, selected by parametric boxes around the corners (addRibCell.m
% does the same). Tip corners at (w/2, +-t/2), root corners at (s/2, +-t/2).
% NOTE addRibCell's tags: fil2 = TIP (r2), fil1 = ROOT (r1).
selTol = '1[nm]';
boxTip = wp.geom.create('fil2_u_box', 'BoxSelection');
boxTip.set('entitydim', 0);
boxTip.set('condition', 'inside');
boxTip.set('xmin', ['w/2-' selTol]);   boxTip.set('xmax', ['w/2+' selTol]);
boxTip.set('ymin', ['-t/2-' selTol]);  boxTip.set('ymax', ['t/2+' selTol]);
filTip = wp.geom.create('fil2_u', 'Fillet');
filTip.set('radius', 'r2');
filTip.selection('point').named('fil2_u_box');

boxRoot = wp.geom.create('fil1_u_box', 'BoxSelection');
boxRoot.set('entitydim', 0);
boxRoot.set('condition', 'inside');
boxRoot.set('xmin', ['s/2-' selTol]);  boxRoot.set('xmax', ['s/2+' selTol]);
boxRoot.set('ymin', ['-t/2-' selTol]); boxRoot.set('ymax', ['t/2+' selTol]);
filRoot = wp.geom.create('fil1_u', 'Fillet');
filRoot.set('radius', 'r1');
filRoot.selection('point').named('fil1_u_box');

if ~isWedge
    beam = geom.create('ribcell', 'Extrude');
    beam.label('Rib unit cell (diamond)');
    beam.selection('input').set({'wp1'});
    beam.setIndex('distance', 'th', 0);
    beam.set('selresult', true);             % -> geom1_ribcell_dom
else
    % The slab as above (its outline carries the fillets), then the wedge
    % under that filleted outline, united into ONE object with the tag
    % 'ribcell', so every name below (and in the drivers) is the same for
    % both cross-sections.
    %
    % CHECK AFTER THE FIRST BUILD (model = build_omc_unitcell_rib('','wedge');
    % mphgeom(model); or mphlaunch(model)):
    %   1. One beam domain; the wedge hangs below the slab (x < -th/2) and
    %      its footprint is exactly the filleted outline (look from -x).
    %   2. Depth below the slab (mphinterp / measure): s/2*tan(30deg) under
    %      the spine centreline between ribs; t/2*tan(30deg) under the rib
    %      centreline mid-arm; wd under the junction centre (z = 0, y = 0),
    %      which is the deepest point; 0 along the whole outline.
    %   3. Work-plane orientation: the rib prism must lie under the rib
    %      (|z| <= t/2), not rotated by 90 deg ('zx' plane: local axes (z, x)).
    %   4. The root cones open DOWNWARDS (radius r1 at the slab, larger
    %      below) and the tip cones point DOWNWARDS (apex below the slab).
    %      If the pointed Cone (rtop = 0) is refused, replace it by a Revolve.
    %   5. The builder's own checks pass (periodic faces pair up, no empty
    %      selection), and mesh.run succeeds.
    slab = geom.create('ribslab', 'Extrude');
    slab.label('Rib slab');
    slab.selection('input').set({'wp1'});
    slab.setIndex('distance', 'th', 0);

    tanA = 'tan(wedge_angle)';
    xBot = '-th/2';                          % bottom face of the slab

    % --- S: spine prism outside the junction zone, |z| >= t/2 + r1 -------
    % Half triangle in the cross-section plane (x, y), as Polygon 1 /
    % Extrude 1 of omc_unitcell_iso.mph with b = s and the y >= 0 half.
    zStart = {'-a/2', 't/2+r1'};
    for i = 1:2
        wps = geom.create(sprintf('wp_ws%d', i), 'WorkPlane');
        wps.set('unite', true);
        wps.set('planetype', 'quick');
        wps.set('quickplane', 'xy');
        wps.set('quickz', zStart{i});
        pol = wps.geom.create('pol', 'Polygon');
        pol.set('source', 'table');
        pol.set('table', {xBot, '0'; ...
                          xBot, 's/2'; ...
                          [xBot '-s/2*' tanA], '0'});
        ws = geom.create(sprintf('wedge_s%d', i), 'Extrude');
        ws.label(sprintf('Spine wedge %d', i));
        ws.selection('input').set({sprintf('wp_ws%d', i)});
        ws.setIndex('distance', 'a/2-t/2-r1', 0);
    end

    % --- A: rib prism, y from s/2 + r1 to the tip -------------------------
    % Triangle in the plane y = s/2 + r1, whose local axes are (z, x) and
    % whose normal is +y. Base: the bottom face of the rib, z = -t/2 .. t/2.
    wpr = geom.create('wp_wr', 'WorkPlane');
    wpr.set('unite', true);
    wpr.set('planetype', 'quick');
    wpr.set('quickplane', 'zx');
    wpr.set('quicky', 's/2+r1');
    pol = wpr.geom.create('pol', 'Polygon');
    pol.set('source', 'table');
    pol.set('table', {'-t/2', xBot; ...
                      't/2',  xBot; ...
                      '0',    [xBot '-t/2*' tanA]});
    wr = geom.create('wedge_r', 'Extrude');
    wr.label('Rib wedge');
    wr.selection('input').set({'wp_wr'});
    wr.setIndex('distance', 'w/2-s/2-r1', 0);

    % --- J: junction zones, one per side of the rib ------------------------
    % Box (0 <= +-z <= t/2+r1, 0 <= y <= s/2+r1, 1.25*wdj deep) minus a cone
    % frustum about the root-arc centre (y = s/2+r1, z = +-(t/2+r1)): radius
    % r1 at the slab, r1 + depth/tan(wedge_angle) below. What is left is
    % depth <= tan(wedge_angle)*(rho - r1).
    hJ = '1.25*wdj';
    zBox = {'0', '-(t/2+r1)'};
    zAxis = {'t/2+r1', '-(t/2+r1)'};
    for i = 1:2
        blk = geom.create(sprintf('jblk%d', i), 'Block');
        blk.set('base', 'corner');
        blk.set('pos',  {[xBot '-' hJ], '0', zBox{i}});
        blk.set('size', {hJ, 's/2+r1', 't/2+r1'});
        con = geom.create(sprintf('jcone%d', i), 'Cone');
        con.set('pos', {[xBot '-' hJ], 's/2+r1', zAxis{i}});
        con.set('axistype', 'x');            % from the deep end up to the slab
        con.set('r', ['r1+' hJ '/' tanA]);   % radius at the deep end
        con.set('h', hJ);
        con.set('specifytop', 'radius');
        con.set('rtop', 'r1');               % radius at the slab = root fillet
        dif = geom.create(sprintf('wedge_j%d', i), 'Difference');
        dif.label(sprintf('Junction wedge %d', i));
        dif.selection('input').set({sprintf('jblk%d', i)});
        dif.selection('input2').set({sprintf('jcone%d', i)});
        dif.set('intbnd', false);
    end

    uni = geom.create('wedge_u', 'Union');
    uni.selection('input').set({'wedge_s1', 'wedge_s2', 'wedge_r', ...
                                'wedge_j1', 'wedge_j2'});
    uni.set('intbnd', false);

    % --- tip plane (hip): keep depth <= tan(wedge_angle)*(w/2 - y) ---------
    wpt = geom.create('wp_wt', 'WorkPlane');
    wpt.set('unite', true);
    wpt.set('planetype', 'quick');
    wpt.set('quickplane', 'xy');
    wpt.set('quickz', '-a/2');
    pol = wpt.geom.create('pol', 'Polygon');
    pol.set('source', 'table');
    pol.set('table', {xBot, '0'; ...
                      xBot, 'w/2'; ...
                      [xBot '-w/2*' tanA], '0'});
    tipb = geom.create('tipblk', 'Extrude');
    tipb.selection('input').set({'wp_wt'});
    tipb.setIndex('distance', 'a', 0);
    hip = geom.create('wedge_h', 'Intersection');
    hip.selection('input').set({'wedge_u', 'tipblk'});
    hip.set('intbnd', false);

    % --- C: tip corners. Remove the corner boxes (they reach into the slab
    % level so no face is coplanar with the wedge), then add quarter cones.
    hC = ['1.5*r2*' tanA];
    zCorner = {'t/2-r2', '-t/2'};            % box start along z
    zCone   = {'t/2-r2', '-(t/2-r2)'};       % tip-arc centres
    for i = 1:2
        for k = 1:2                          % one box to cut, one to intersect
            blk = geom.create(sprintf('cblk%d%d', i, k), 'Block');
            blk.set('base', 'corner');
            blk.set('pos',  {[xBot '-' hC], 'w/2-r2', zCorner{i}});
            blk.set('size', {[hC '+th/4'], 'r2', 'r2'});
        end
        con = geom.create(sprintf('ccone%d', i), 'Cone');
        con.set('pos', {xBot, 'w/2-r2', zCone{i}});
        con.set('axistype', 'cartesian');
        con.set('axis', {'-1', '0', '0'});   % pointing down, away from the slab
        con.set('r', 'r2');                  % radius at the slab = tip fillet
        con.set('h', ['r2*' tanA]);
        con.set('specifytop', 'radius');
        con.set('rtop', '0');                % a point: the apex
        qc = geom.create(sprintf('wedge_c%d', i), 'Intersection');
        qc.label(sprintf('Tip corner wedge %d', i));
        qc.selection('input').set({sprintf('ccone%d', i), sprintf('cblk%d2', i)});
        qc.set('intbnd', false);
    end
    cut = geom.create('wedge_hc', 'Difference');
    cut.selection('input').set({'wedge_h'});
    cut.selection('input2').set({'cblk11', 'cblk21'});
    cut.set('intbnd', false);

    beam = geom.create('ribcell', 'Union');
    beam.label('Rib unit cell with wedge (diamond)');
    beam.selection('input').set({'ribslab', 'wedge_hc', 'wedge_c1', 'wedge_c2'});
    beam.set('intbnd', false);
    beam.set('selresult', true);             % -> geom1_ribcell_dom
end

% Air half-disk: a cylinder of radius air_r along the beam (axis z through
% x = air_xc, y = 0), one period long, with its y < 0 half removed (as
% Block 3 / Difference 2 of omc_unitcell.mph remove that half), minus the
% beam, so air is its own object. Shape and role are those of Work Plane 2 /
% Revolve 1 of omc_unitcell_iso.mph. A Cylinder is used instead of a Revolve:
% its end faces are one face each and it has no edge along the axis, which a
% Revolve would leave inside (or, for a deep wedge, just under) the beam.
% The side of a Cylinder comes as four faces with seams at 0, 90, 180 and
% 270 deg, so the y >= 0 half keeps two curved faces.
if isWedge
    xBeamLo = '-th/2-wd';
else
    xBeamLo = '-th/2';
end
airCyl = geom.create('aircyl', 'Cylinder');
airCyl.set('axistype', 'z');
airCyl.set('pos', {'air_xc', '0', '-a/2'});
airCyl.set('r', 'air_r');
airCyl.set('h', 'a');
airCut = geom.create('aircut', 'Block');     % the y < 0 half, oversized
airCut.set('base', 'corner');
airCut.set('pos',  {'air_xc-1.5*air_r', '-1.5*air_r', '-a'});
airCut.set('size', {'3*air_r', '1.5*air_r', '2*a'});
airHalf = geom.create('airdsk', 'Difference');
airHalf.label('Air half-disk');
airHalf.selection('input').set({'aircyl'});
airHalf.selection('input2').set({'aircut'});
airHalf.set('intbnd', false);

air = geom.create('airdif', 'Difference');
air.label('Air');
air.selection('input').set({'airdsk'});
air.selection('input2').set({'ribcell'});
air.set('keepsubtract', true);
air.set('selresult', true);                  % -> geom1_airdif_dom

% Boundary selections. Bounds are parameter expressions, so they follow the
% geometry when a parameter changes. Unset bounds stay at -Inf / +Inf.
tol = '1[nm]';
xBeam = {'xmin', [xBeamLo '-' tol], 'xmax', ['th/2+' tol]};
yBeam = {'ymin', ['-' tol],      'ymax', ['w/2+' tol]};
zLo   = {'zmin', ['-a/2-' tol],  'zmax', ['-a/2+' tol]};
zHi   = {'zmin', ['a/2-' tol],   'zmax', ['a/2+' tol]};
yZero = {'ymin', ['-' tol],      'ymax', tol};
boxes = { ...
    'zloS',         [xBeam, yBeam, zLo]; ...   % beam face z = -a/2
    'zhiS',         [xBeam, yBeam, zHi]; ...   % beam face z = +a/2
    'zloA',         zLo; ...                   % all faces z = -a/2
    'zhiA',         zHi; ...                   % all faces z = +a/2
    'selYsymSolid', [xBeam, yZero]; ...        % beam face y = 0
    'selYsymAll',   yZero; ...                 % all faces y = 0
    'allA',         {}; ...                    % every face (no bounds)
    'beamA',        [xBeam, yBeam]};           % every face of the beam
for i = 1:size(boxes, 1)
    box = geom.create(boxes{i, 1}, 'BoxSelection');
    box.set('entitydim', 2);
    box.set('condition', 'allvertices');
    bounds = boxes{i, 2};
    for j = 1:2:numel(bounds)
        box.set(bounds{j}, bounds{j + 1});
    end
end
% Probes for the check of selScat below: boxes on the curved air surface at
% 45 and 135 deg about the axis, mid-cell, which only that surface reaches
% (as airOuterBox1/2 of omc-comsol-chang/buildRibUnitCell_LN.m). They need
% air_gap > 0.08*air_r, which the defaults give with a wide margin.
probeX = {'air_xc+air_r/sqrt(2)', 'air_xc-air_r/sqrt(2)'};
for i = 1:2
    box = geom.create(sprintf('scatP%d', i), 'BoxSelection');
    box.set('entitydim', 2);
    box.set('condition', 'intersects');
    box.set('xmin', [probeX{i} '-0.05*air_r']);
    box.set('xmax', [probeX{i} '+0.05*air_r']);
    box.set('ymin', 'air_r/sqrt(2)-0.05*air_r');
    box.set('ymax', 'air_r/sqrt(2)+0.05*air_r');
    box.set('zmin', '-a/4');
    box.set('zmax', 'a/4');
end
unions = { ...
    'selPerSolid', {'zloS', 'zhiS'}; ...           % Floquet pair, beam
    'selPerAll',   {'zloA', 'zhiA'}; ...           % Floquet pair, beam + air
    'scatProbe',   {'scatP1', 'scatP2'}};          % check of selScat only
for i = 1:size(unions, 1)
    sel = geom.create(unions{i, 1}, 'UnionSelection');
    sel.set('entitydim', 2);
    sel.set('input', unions{i, 2});
end
% Outer air surface (scattering): every face that is not on z = +-a/2, not
% on y = 0 and not a face of the beam. What is left is the curved surface of
% the half-disk, whatever its radius, so the selection follows air_r.
scat = geom.create('selScat', 'DifferenceSelection');
scat.set('entitydim', 2);
scat.set('add', {'allA'});
scat.set('subtract', {'zloA', 'zhiA', 'selYsymAll', 'beamA'});

geom.run;                                    % run (not runAll) finalizes

% An empty selection would give a model that solves but has a missing
% boundary condition, so stop here instead.
needed = {'geom1_ribcell_dom', 'geom1_airdif_dom', 'geom1_selPerSolid', ...
          'geom1_selPerAll', 'geom1_selYsymSolid', 'geom1_selYsymAll', ...
          'geom1_selScat'};
for i = 1:numel(needed)
    if isempty(model.selection(needed{i}).entities())
        error('build_omc_unitcell_rib:emptySelection', ...
            'Selection %s is empty for the default parameters.', needed{i});
    end
end
nPer = numel(model.selection('geom1_selPerSolid').entities());
if ~isWedge && nPer ~= 2
    error('build_omc_unitcell_rib:periodicFaces', ...
        'Expected 2 periodic beam faces (z = -a/2 and z = +a/2).');
end
% With the wedge each end of the beam may be one face or two (slab and
% wedge part): the same number is needed at z = -a/2 and at z = +a/2.
nLo = numel(model.selection('geom1_zloS').entities());
nHi = numel(model.selection('geom1_zhiS').entities());
if isWedge && (nLo < 1 || nLo ~= nHi || nPer ~= nLo + nHi)
    error('build_omc_unitcell_rib:periodicFaces', ...
        'Periodic beam faces do not pair up: %d at z = -a/2, %d at z = +a/2.', ...
        nLo, nHi);
end
% Beam + air end faces (ewfd pc1): the same number at both ends, and more
% than the beam's own (one half-disk face of air per end is expected).
nLoA = numel(model.selection('geom1_zloA').entities());
nHiA = numel(model.selection('geom1_zhiA').entities());
if nLoA ~= nHiA || nLoA <= nLo
    error('build_omc_unitcell_rib:periodicFacesAir', ...
        ['Periodic beam + air faces do not pair up: %d at z = -a/2 (%d of ' ...
         'the beam), %d at z = +a/2 (%d of the beam).'], nLoA, nLo, nHiA, nHi);
end
% selScat must be exactly the curved air surface: the faces found by the two
% probes on it. A face missing from selScat would keep the default PEC.
scatFaces  = sort(double(model.selection('geom1_selScat').entities()));
probeFaces = sort(double(model.selection('geom1_scatProbe').entities()));
if ~isequal(scatFaces(:), probeFaces(:))
    error('build_omc_unitcell_rib:scatFaces', ...
        ['selScat holds faces %s, but the curved air surface is faces %s ' ...
         '(probes scatP1, scatP2).'], ...
        mat2str(scatFaces(:).'), mat2str(probeFaces(:).'));
end
if numel(scatFaces) ~= 2
    warning('build_omc_unitcell_rib:scatFaceCount', ...
        'Expected 2 curved air faces (two quarter cylinders), found %d.', ...
        numel(scatFaces));
end

%% Definitions: integration over the beam (mode classification)
intop = comp.cpl.create('intop1', 'Integration');
intop.selection.named('geom1_ribcell_dom');

%% Materials (values copied from omc_unitcell.mph, except the elasticity)
% ANISOTROPIC (cubic) diamond, unlike omc_unitcell.mph / omc_unitcell_iso.mph,
% which use isotropic E = 1050 GPa, nu = 0.1 (the <100> values). Stiffness
% C11/C12/C44 are global parameters (configs/materials.yaml). Voigt order
% (xx, yy, zz, yz, xz, xy) and DVo layout (upper triangle, column by column)
% as omc-comsol-chang/LoadMaterialParams.m / runBands.m.
%
% Crystal rotation rxt, as omc-comsol-chang/RotateXtalTensor.m (P.rxtal):
% there the frame is beam x, width y, thickness z and the crystal is rotated
% about z (the slab normal). Here the frame is beam z, width y, thickness x,
% i.e. the Voigt permutation [3 2 1 6 5 4] of (xx, yy, zz, yz, xz, xy), and
% the rotation is about x. For a cubic crystal, with Cani = C11-C12-2*C44,
% dCrot = Cani/2*sin(2*rxt)^2 and gCrot = Cani/4*sin(4*rxt):
%   D22 = D33 = C11-dCrot,  D23 = C12+dCrot,  D44 = C44+dCrot,
%   D24 = gCrot,  D34 = -gCrot,  all others as the unrotated cube axes.
% rxt = 0: the model axes are the cube axes, beam along <100>. rxt = 45[deg]:
% beam and width along <110>. gCrot = 0 at multiples of 45 deg only; any other
% rxt would couple y-even and y-odd strain and break the y = 0 mirror of the
% half model, so only 0 and 45[deg] are used (python-scripts resolve_rxtal).
% The label must keep starting with "Diamond": set_diamond_index and
% patch_hole_templates_index.m find this material by that prefix.
diamond = comp.material.create('mat1', 'Common');
diamond.label('Diamond (cubic, rotated by rxt about x)');
diamond.selection.named('geom1_ribcell_dom');
def = diamond.propertyGroup('def');
def.set('thermalexpansioncoefficient', {'0.8e-6[1/K]'});
def.set('heatcapacity', '520[J/(kg*K)]');
def.set('relpermittivity', {'n_diamond^2'});   % n = 2.386 (was 5.76 = 2.40^2)
def.set('density', '3515[kg/m^3]');
def.set('thermalconductivity', {'990[W/(m*K)]'});
aniso = diamond.propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
aniso.set('DVo', {'C11', ...
                  'C12', 'C11-dCrot', ...
                  'C12', 'C12+dCrot', 'C11-dCrot', ...
                  '0',   'gCrot',     '-gCrot',    'C44+dCrot', ...
                  '0',   '0',         '0',         '0',   'C44', ...
                  '0',   '0',         '0',         '0',   '0',   'C44'});

vacuum = comp.material.create('mat2', 'Common');
vacuum.label('Perfect vacuum');
vacuum.selection.named('geom1_airdif_dom');
vdef = vacuum.propertyGroup('def');
vdef.set('density', '0[kg/m^3]');
vdef.set('relpermeability', {'1'});
vdef.set('relpermittivity', {'1'});
vdef.set('electricconductivity', {'0[S/m]'});

%% Solid Mechanics (beam only; anisotropic linear elastic, D from mat1)
% D comes from mat1, already rotated by rxt (above), so lemm1 keeps the
% Global coordinate system: no rotated system is needed here.
solid = comp.physics.create('solid', 'SolidMechanics', 'geom1');
solid.selection.named('geom1_ribcell_dom');
lemm = solid.feature('lemm1');
lemm.set('SolidModel', 'Anisotropic');
lemm.set('AnisotropicOption', 'AnisotropicVo');
lemm.set('DVo_mat', 'from_mat');

pc = solid.create('pc1', 'PeriodicCondition', 2);
pc.selection.named('geom1_selPerSolid');
pc.set('PeriodicType', 'Floquet');
pc.set('kFloquet', {'0'; '0'; 'kF'});

sym = solid.create('sym1', 'SymmetrySolid', 2);       % used by std1
sym.selection.named('geom1_selYsymSolid');
asym = solid.create('as1', 'Antisymmetry', 2);        % used by std3
asym.selection.named('geom1_selYsymSolid');

%% Electromagnetic Waves, Frequency Domain (beam + air)
ewfd = comp.physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');

epc = ewfd.create('pc1', 'PeriodicCondition', 2);
epc.selection.named('geom1_selPerAll');
epc.set('PeriodicType', 'Floquet');
epc.set('Floquet_source', 'UserDefined');
epc.set('kFloquet', {'0'; '0'; 'kF'});

sctr = ewfd.create('sctr1', 'Scattering', 2);
sctr.selection.named('geom1_selScat');
pmc = ewfd.create('pmc1', 'PerfectMagneticConductor', 2);   % used by std4
pmc.selection.named('geom1_selYsymAll');
pec = ewfd.create('pec2', 'PerfectElectricConductor', 2);   % used by std2
pec.selection.named('geom1_selYsymAll');

%% Mesh: physics-controlled, element size Normal (as in the template)
% For ewfd this sets the maximum element size from the wavelength of the
% study's shift: (1550 nm)/5 = 310 nm in air and 310 nm/n in the beam (read
% from the block-air omc_unitcell_rib.mph: Size 3.1e-7, Size 1 1.299e-7), so
% the air half-disk gets no size node of its own.
% TODO(confirm): mesh convergence of the optical bands; nothing changed here.
mesh = comp.mesh.create('mesh1');
mesh.autoMeshSize(5);

%% Studies: tag, label, active physics, feature disabled in that study
studies = { ...
    'std1', 'mech sym',     'solid', 'as1'; ...
    'std3', 'mech antisym', 'solid', 'sym1'; ...
    'std2', 'opt TE',       'ewfd',  'pmc1'; ...
    'std4', 'opt TM',       'ewfd',  'pec2'};
for i = 1:size(studies, 1)
    study = model.study.create(studies{i, 1});
    study.label(studies{i, 2});
    eig = study.create('eig', 'Eigenfrequency');
    if strcmp(studies{i, 3}, 'solid')
        eig.set('neigsactive', true);
        eig.set('neigs', nEigMech);
        eig.set('eigunit', 'GHz');
        eig.set('shiftactive', true);
        eig.set('shift', '8');
        eig.activate('ewfd', false);
    else
        eig.set('eigunit', 'THz');
        eig.set('shiftactive', true);
        eig.set('shift', 'c_const/1.55[um]');
        eig.activate('solid', false);
    end
    eig.set('useadvanceddisable', true);
    eig.set('disabledphysics', {[studies{i, 3} '/' studies{i, 4}]});
end

% Create the solver sequences now, in this order, so each study's dataset is
% fixed: dset1 = mech sym, dset2 = mech antisym, dset3 = opt TE, dset4 = opt TM.
for i = 1:size(studies, 1)
    model.study(studies{i, 1}).createAutoSequences('sol');
end

mesh.run;

%% Save
if ~isempty(outFile)
    mphsave(model, outFile);
    fprintf('Saved %s\n', outFile);
end
end
