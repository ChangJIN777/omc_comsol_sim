function model = build_omc_unitcell_rib(outFile)
% BUILD_OMC_UNITCELL_RIB  Build the parametric rib unit-cell template
% (omc_unitcell_rib.mph) via COMSOL LiveLink for MATLAB.
%
%   model = build_omc_unitcell_rib()         % saves comsol/omc_unitcell_rib.mph
%   model = build_omc_unitcell_rib(outFile)  % saves to outFile
%   model = build_omc_unitcell_rib('')       % build in memory only, no save
%
% Rib counterpart of omc_unitcell.mph. The node structure (tags, labels,
% physics, boundary conditions, material, mesh approach, studies) is copied
% from omc_unitcell.mph; the differences are the geometry, 20 (not 12)
% eigenfrequencies in the two mechanical studies, and ANISOTROPIC (cubic)
% diamond with C11/C12/C44 as global parameters, beam along <100>. The geometry is the
% recipe of omc-comsol-chang/buildRibUnitCell.m (work plane -> half-width
% rectangle -> polygon void -> Compose -> Extrude), with two changes:
%   * every dimension is a parameter EXPRESSION, so the saved .mph is tunable
%   * the frame is the template's: beam axis z, width y, thickness x
%     (buildRibUnitCell uses x, y, z), so Floquet k = (0,0,kF) is unchanged
%
% Global parameters (all referenced by name in the geometry):
%   a   period along z                      th  slab thickness (x)
%   w   full beam width (y); the model      ai  spine width: the void starts
%       draws the y >= 0 half                   at y = ai/2
%   wo  void length along z                 ho  void height along y
%   wi  tooth length along z (wi < wo)      hi  tooth height along y (hi < ho)
%   kF  Floquet wavevector along z (swept by the driver)
%
% Half cell seen from +x (y up, z right):
%
%   y = w/2      +---------------------------+
%                |        diamond            |
%   ai/2 + ho    |    +-----------------+    |
%                |    |      void       |    |
%   ai/2 + hi    |    |   +---------+   |    |
%                |    |   |  tooth  |   |    |
%   ai/2         |    +---+         +---+    |
%                |          spine            |
%   y = 0        +---------------------------+   <- y symmetry plane
%              -a/2  -wo/2 -wi/2  wi/2 wo/2  a/2     (z)
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
% Deliberate differences from buildRibUnitCell.m:
%   * one period only. P.nperiod / P.holeatedge are dropped: the builder's
%     hole loop never offsets its polygons, and a Floquet cell is one period.
%   * no thickness-symmetry cut (P.mbevenz). omc_unitcell.mph has none.
%   * debug figure / mphgeom / mphplot / disp calls removed.
%   * selections are parametric boxes that follow the drawn cell
%     (z = -a/2..a/2, y = 0..w/2), and the geometry is finalized with run.
%
% Requires: COMSOL 6.2+ with LiveLink for MATLAB connected (mphstart), the
% Structural Mechanics (or MEMS / Acoustics) and Wave Optics modules.
%
% See also OMC_MECHANICAL_LIVELINK_RIB, OMC_MECHANICAL_LIVELINK.

import com.comsol.model.*
import com.comsol.model.util.*

if nargin < 1
    outFile = fullfile(fileparts(mfilename('fullpath')), 'omc_unitcell_rib.mph');
end
if ~isempty(outFile) && exist(outFile, 'file')
    error('build_omc_unitcell_rib:fileExists', ...
        '%s already exists. Delete it or pass another file name.', outFile);
end

nEigMech = 20;      % eigenfrequencies per mechanical study (template: 12)

%% Model and global parameters
model = ModelUtil.create(char(ModelUtil.uniquetag('Model')));
model.label('omc_unitcell_rib.mph');

% Defaults: omc-comsol-chang/test_nanobeamRibGeometry.m ('2D_ribs').
% TODO(confirm) values.
params = { ...
    'a',  '1500[nm]', 'lattice constant (period along z)'; ...
    'w',  '3000[nm]', 'full beam width (y); the y >= 0 half is drawn'; ...
    'th', '160[nm]',  'beam thickness (x)'; ...
    'wi', '670[nm]',  'tooth (inner) length along z'; ...
    'wo', '1000[nm]', 'void (outer) length along z'; ...
    'hi', '955[nm]',  'tooth (inner) height along y'; ...
    'ho', '975[nm]',  'void (outer) height along y'; ...
    'ai', '150[nm]',  'spine width; the void starts at y = ai/2'; ...
    'kF', 'pi/a',     'Floquet wavevector along z (swept by driver)'; ...
    'C11', '1076[GPa]', 'diamond cubic stiffness C11 (configs/materials.yaml)'; ...
    'C12', '125[GPa]',  'diamond cubic stiffness C12 (configs/materials.yaml)'; ...
    'C44', '577[GPa]',  'diamond cubic stiffness C44 (configs/materials.yaml)'};
for i = 1:size(params, 1)
    model.param.set(params{i, 1}, params{i, 2}, params{i, 3});
end

%% Geometry (buildRibUnitCell recipe, template frame)
comp = model.component.create('comp1', true);
geom = comp.geom.create('geom1', 3);
geom.lengthUnit('m');

% Work plane x = -th/2. Its local axes are (y, z); the normal is +x.
wp = geom.create('wp1', 'WorkPlane');
wp.set('unite', true);
wp.set('planetype', 'quick');
wp.set('quickplane', 'yz');
wp.set('quickx', '-th/2');

% y-symmetric half of the slab: y = 0..w/2, z = -a/2..a/2
slab = wp.geom.create('ucellplane', 'Rectangle');
slab.set('type', 'solid');
slab.set('base', 'corner');
slab.set('pos',  {'0', '-a/2'});
slab.set('size', {'w/2', 'a'});

% Void polygon. Same vertices and order as buildRibUnitCell.m, written as
% (y, z) expressions. Its repeated closing vertex is left out.
void = wp.geom.create('hole_1_upper', 'Polygon');
void.set('source', 'table');
void.set('table', { ...
    'ai/2',    '-wo/2'; ...
    'ai/2',    '-wi/2'; ...
    'ai/2+hi', '-wi/2'; ...
    'ai/2+hi', 'wi/2'; ...
    'ai/2',    'wi/2'; ...
    'ai/2',    'wo/2'; ...
    'ai/2+ho', 'wo/2'; ...
    'ai/2+ho', '-wo/2'});

cut = wp.geom.create('ucellComp', 'Compose');
cut.selection('input').set({'ucellplane', 'hole_1_upper'});
cut.set('formula', 'ucellplane-hole_1_upper');

beam = geom.create('ucellholes', 'Extrude');
beam.label('Rib unit cell (diamond)');
beam.selection('input').set({'wp1'});
beam.setIndex('distance', 'th', 0);
beam.set('selresult', true);                 % -> geom1_ucellholes_dom

% Air box with the template's extents (5*th in x, 1.5*w in y after the
% y < 0 half is removed), minus the beam, so air is its own object.
airBox = geom.create('airblk', 'Block');
airBox.set('base', 'corner');
airBox.set('pos',  {'-2.5*th', '0', '-a/2'});
airBox.set('size', {'5*th', '1.5*w', 'a'});

air = geom.create('airdif', 'Difference');
air.label('Air');
air.selection('input').set({'airblk'});
air.selection('input2').set({'ucellholes'});
air.set('keepsubtract', true);
air.set('selresult', true);                  % -> geom1_airdif_dom

% Boundary selections. Bounds are parameter expressions, so they follow the
% geometry when a parameter changes. Unset bounds stay at -Inf / +Inf.
tol = '1[nm]';
xBeam = {'xmin', ['-th/2-' tol], 'xmax', ['th/2+' tol]};
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
    'xloA', {'xmin', ['-2.5*th-' tol], 'xmax', ['-2.5*th+' tol]}; ...
    'xhiA', {'xmin', ['2.5*th-' tol],  'xmax', ['2.5*th+' tol]}; ...
    'yhiA', {'ymin', ['1.5*w-' tol],   'ymax', ['1.5*w+' tol]}};
for i = 1:size(boxes, 1)
    box = geom.create(boxes{i, 1}, 'BoxSelection');
    box.set('entitydim', 2);
    box.set('condition', 'allvertices');
    bounds = boxes{i, 2};
    for j = 1:2:numel(bounds)
        box.set(bounds{j}, bounds{j + 1});
    end
end
unions = { ...
    'selPerSolid', {'zloS', 'zhiS'}; ...           % Floquet pair, beam
    'selPerAll',   {'zloA', 'zhiA'}; ...           % Floquet pair, beam + air
    'selScat',     {'xloA', 'xhiA', 'yhiA'}};      % outer air faces
for i = 1:size(unions, 1)
    sel = geom.create(unions{i, 1}, 'UnionSelection');
    sel.set('entitydim', 2);
    sel.set('input', unions{i, 2});
end

geom.run;                                    % run (not runAll) finalizes

% An empty selection would give a model that solves but has a missing
% boundary condition, so stop here instead.
needed = {'geom1_ucellholes_dom', 'geom1_airdif_dom', 'geom1_selPerSolid', ...
          'geom1_selPerAll', 'geom1_selYsymSolid', 'geom1_selYsymAll', ...
          'geom1_selScat'};
for i = 1:numel(needed)
    if isempty(model.selection(needed{i}).entities())
        error('build_omc_unitcell_rib:emptySelection', ...
            'Selection %s is empty for the default parameters.', needed{i});
    end
end
if numel(model.selection('geom1_selPerSolid').entities()) ~= 2
    error('build_omc_unitcell_rib:periodicFaces', ...
        'Expected 2 periodic beam faces (z = -a/2 and z = +a/2).');
end

%% Definitions: integration over the beam (mode classification)
intop = comp.cpl.create('intop1', 'Integration');
intop.selection.named('geom1_ucellholes_dom');

%% Materials (values copied from omc_unitcell.mph, except the elasticity)
% ANISOTROPIC (cubic) diamond, unlike omc_unitcell.mph / omc_unitcell_iso.mph,
% which use isotropic E = 1050 GPa, nu = 0.1 (the <100> values). Stiffness
% C11/C12/C44 are global parameters (configs/materials.yaml). The model axes
% (thickness x, width y, beam z) are the crystal cube axes, so the tensor is
% unrotated: beam along <100>. Voigt order and DVo layout as
% omc-comsol-chang/LoadMaterialParams.m / runBands.m.
diamond = comp.material.create('mat1', 'Common');
diamond.label('Diamond (100)');
diamond.selection.named('geom1_ucellholes_dom');
def = diamond.propertyGroup('def');
def.set('thermalexpansioncoefficient', {'0.8e-6[1/K]'});
def.set('heatcapacity', '520[J/(kg*K)]');
def.set('relpermittivity', {'5.76'});
def.set('density', '3515[kg/m^3]');
def.set('thermalconductivity', {'990[W/(m*K)]'});
aniso = diamond.propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
aniso.set('DVo', {'C11', ...
                  'C12', 'C11', ...
                  'C12', 'C12', 'C11', ...
                  '0',   '0',   '0',   'C44', ...
                  '0',   '0',   '0',   '0',   'C44', ...
                  '0',   '0',   '0',   '0',   '0',   'C44'});

vacuum = comp.material.create('mat2', 'Common');
vacuum.label('Perfect vacuum');
vacuum.selection.named('geom1_airdif_dom');
vdef = vacuum.propertyGroup('def');
vdef.set('density', '0[kg/m^3]');
vdef.set('relpermeability', {'1'});
vdef.set('relpermittivity', {'1'});
vdef.set('electricconductivity', {'0[S/m]'});

%% Solid Mechanics (beam only; anisotropic linear elastic, D from mat1)
solid = comp.physics.create('solid', 'SolidMechanics', 'geom1');
solid.selection.named('geom1_ucellholes_dom');
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
