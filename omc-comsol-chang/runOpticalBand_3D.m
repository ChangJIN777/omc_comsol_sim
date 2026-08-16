function ds = runOpticalBand_3D(P)

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

ModelUtil.clear();
clear model     
model = ModelUtil.create('model');

beamMat = P.beamMat;
a = P.a;
th = P.th;

kpts = P.kpts;
nbands = P.nbands;
freq = P.optical_freq;
max_dof = P.max_dof;
meshSize = P.meshSize;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric


% create P.datLoc if it does not already exist -- normalise separators first.
% The test scripts hardcode '\', which on macOS/Linux is an ordinary filename
% character rather than a separator, so without this mkdir makes one oddly-named
% folder instead of the intended tree. Mirrors what runBands_2D does. Needed
% HERE, before the plotgeom block below, because that block is the first thing
% in this function to write a file - the saveplots mkdir near the end is far too
% late, and is only for the per-design subfolder.
if isfield(P,'datLoc') && ~isempty(P.datLoc)
    P.datLoc = strrep(strrep(P.datLoc,'\',filesep),'/',filesep);
    if ~strcmp(P.datLoc(end),filesep)
        P.datLoc = [P.datLoc,filesep];
    end
    if ~exist(P.datLoc,'dir')
        mkdir(P.datLoc)
    end
end

%% Define k-points for sweep over wavevectors (2D band structure)
% adapted from phononic crystal model on COMSOL

% parameter node for COMSOL model
% k runs from 0 to 3: 0-->1 for Gamma-X, 1-->2 for X-->M, 2-->3 for
% M-Gamma
model.param.set('k', '0');
model.param.set('a', [num2str(a),'[m]']);

if strcmp(P.unitcell,'hexagonal')
    model.param.set('kx', 'if(k<1,(pi/a)*k*(sqrt(3)/2),if(k<2,(sqrt(3)/2)*pi/a,(sqrt(3)/2)*(3-k)*pi/a))');
    model.param.set('ky', 'if(k<1,(pi/a)*k*(-1/2),if(k<2,(k-1-1/2)*pi/a,(1/2)*(3-k)*pi/a))');
    
    for ki = 0:3*kpts-1
        ds.k_norm(ki+1,1) = ki/kpts;
        ds.kx_norm(ki+1,1) = (((sqrt(3)/2)*ki/kpts)*(ki<kpts)+...                  % Gamma-X
                            (sqrt(3)/2)*(ki>=kpts && ki<2*kpts)+...              % X-M
                            (sqrt(3)/2)*(3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
        % The Gamma-X term carries ki/kpts because COMSOL's ky above ramps as
        % (pi/a)*k*(-1/2) along this leg. A bare -1/2 recorded ky_norm = -1/2 at
        % the Gamma point, where the solve actually used ky = 0. Harmless while
        % only kx_norm fed the light line in solveOpticalBands; wrong as soon as
        % ky_norm does.
        ds.ky_norm(ki+1,1) = ((-1/2)*(ki/kpts)*(ki<kpts)+...                % Gamma-X
                            ((ki-kpts)/kpts-1/2)*(ki>=kpts && ki<2*kpts)+... % X-M
                            (1/2)*(3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    end
    
    % compile expressions for input to COMSOL model
    kliststr = ['range(0,1/',num2str(kpts),',3-1/',num2str(kpts),')'];
    for ki = 1:3*kpts
        kparamstr{ki} = ['"k", "',num2str((ki-1)/kpts),'"'];
    end
else
    model.param.set('kx', 'if(k<1,pi/a*k,if(k<2,pi/a,(3-k)*pi/a))');
    model.param.set('ky', 'if(k<1,0,if(k<2,(k-1)*pi/a,(3-k)*pi/a))');
    
    for ki = 0:3*kpts-1
        ds.k_norm(ki+1,1) = ki/kpts;
        ds.kx_norm(ki+1,1) = ((ki/kpts)*(ki<kpts)+...                  % Gamma-X
                            1*(ki>=kpts && ki<2*kpts)+...              % X-M
                            (3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
        ds.ky_norm(ki+1,1) = (0*(ki<kpts)+...                          % Gamma-X
                            (ki-kpts)/kpts*(ki>=kpts && ki<2*kpts)+... % X-M
                            (3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    end
    % compile expressions for input to COMSOL model
    kliststr = ['range(0,1/',num2str(kpts),',3-1/',num2str(kpts),')'];
    for ki = 1:3*kpts
        kparamstr{ki} = ['"k", "',num2str((ki-1)/kpts),'"'];
    end
end 

%% Set up the geometry
if strcmp(P.celltype,'cross')
    [model,P] = DrawCrossUnitCell(model,P);
elseif strcmp(P.celltype,'boomerang')
    [model,P] = buildBoomerangUnitCell(model,P);
elseif strcmp(P.celltype,'boomerang_lower')
    [model,P] = buildLowerBoomerangUnitCell_2D(model,P);
elseif strcmp(P.celltype,'snowflake')
    [model,P] = buildSnowflakeUnitCell_optical(model,P);
elseif strcmp(P.celltype,'hole')
    [model,P] = buildHoleUnitCell_2D(model,P);
elseif strcmp(P.celltype,'hole_strip')
    [model,P] = buildHoleStrip_3D(model,P);
elseif strcmp(P.celltype,'rib')
    [model,P] = buildRibUnitCell_LN(model,P);
end

if P.plotgeom
    figure;
    mphgeom(model)
    % P.datLoc already carries a trailing separator from the block at the top,
    % so no separator is inserted here. The previous [P.datLoc,'\',...] both
    % doubled it and hardcoded a backslash.
    pathFig = [P.datLoc,P.fileBase,'_geom'];
    saveas(gcf,[pathFig,'.fig']);
    saveas(gcf,[pathFig,'.png']);
end

%% Define material and properties
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat1').selection.set([1 2]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').selection.set([2]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
if strcmp(P.celltype,'rib')
    model.component('comp1').material('mat2').selection.set([1]);
    model.component('comp1').material('mat1').selection.set([2 3 4 5 6 7 8 9 10]);
end

model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'2.406' '0' '0' '0' '2.406' '0' '0' '0' '2.406'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
% Air: zero conductivity, for the reason spelled out at the second (winning)
% assignment to this same property below.
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});

if strcmp(P.beamMat,'diamond')
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.406' '0' '0' '0' '2.406' '0' '0' '0' '2.406'});
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
elseif strcmp(P.beamMat,'diamond_telecom')
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.386' '0' '0' '0' '2.386' '0' '0' '0' '2.386'});
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
elseif strcmp(P.beamMat,'SiC')
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.5' '0' '0' '0' '2.5' '0' '0' '0' '2.5'});
elseif strcmp(P.beamMat,'LN')
    model.component('comp1').material('mat1').label('LiNbO3 (Lithium niobate) (Zelmon et al. 1997: n(o) 0.4-5.0 um)');
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', '');
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', '');
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func.create('an1', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func.create('an2', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').label('Sellmeyer formula - extraordinary refractive index');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').set('funcname', 'neref_sel');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').set('expr', '(1 + ((lambda0^(2))*2.9804e12)/((lambda0^(2))*1e12 - 0.02047) + ((lambda0^(2))*0.5981e12)/((lambda0^(2))*1e12 - 0.0666) + ((lambda0^(2))*8.9543e12)/((lambda0^(2))*1e12 - 416.08))^(1./2.)');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').set('args', {'lambda0'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').set('argunit', {'m'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an1').set('plotargs', {'lambda0' '0.4[um]' '5[um]'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').label('Sellmeyer formula - ordinary refractive index');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').set('funcname', 'noref_sel');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').set('expr', '(1 + ((lambda0^(2))*2.6734e12)/((lambda0^(2))*1e12 - 0.01764) + ((lambda0^(2))*1.22901e12)/((lambda0^(2))*1e12 - 0.05914) + ((lambda0^(2))*12.614e12)/((lambda0^(2))*1e12 - 474.60))^(1./2.)');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').set('args', {'lambda0'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').set('argunit', {'m'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').func('an2').set('plotargs', {'lambda0' '0.4[um]' '5[um]'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'neref_sel(lbd0)' '0' '0' '0' 'noref_sel(lbd0)' '0' '0' '0' 'noref_sel(lbd0)'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').addInput('frequency');
else
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'3.5' '0' '0' '0' '3.5' '0' '0' '0' '3.5'});
end
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
% Exactly zero, not a small regularizing value. mat2 is air, which carries the
% Scattering condition (sctr1) below, and that condition's kinp contains
% sigma/(epsr*iomega*eps0). At sigma = 0 COMSOL cancels the term symbolically,
% so it cannot go singular however the eigensolver linearizes; at 1e-12 the
% term survives and divides by zero whenever the linearization point is 0.
% Belt and braces with eigref on the eigenvalue solver further down.
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});

%% setup the physics and the boundary conditions
% Out-of-plane radiation is absorbed by a low-reflecting SCATTERING boundary on
% the top of the air region, rather than by a PML. The PML that used to sit here
% is kept commented rather than deleted so it is easy to put back:
%
% model.component('comp1').coordSystem.create('pml1', 'PML');
% model.component('comp1').coordSystem('pml1').selection.set([3]);
%
% Why the swap is more than cosmetic: the PML was pinned to the hard-coded domain
% index 3, which only means "the absorbing slab" for whichever builder happened
% to produce the domain ordering it was written against. The scattering condition
% instead selects its boundary geometrically - see airTopBoundary below - so it
% follows the geometry rather than an index. The trade-off is accuracy: a
% first-order low-reflecting condition is only near-perfect for normal
% incidence, so it reflects more than a well-tuned PML at grazing angles. Put the
% air region a few wavelengths above the slab to keep that reflection off the
% guided modes, and treat Q factors from this model as a lower bound.
model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('emw').feature('pc1').selection.named('geom1_xboundaries_bnd');
model.component('comp1').physics('emw').create('pc2', 'PeriodicCondition', 2);
model.component('comp1').physics('emw').feature('pc2').selection.named('geom1_yboundaries_bnd');

% Low-reflecting boundary on the air region's top face. The boundary is located
% by normal direction and z offset rather than by a literal index.
sctrBnd = airTopBoundary(model,P);
model.component('comp1').physics('emw').create('sctr1', 'Scattering', 2);
model.component('comp1').physics('emw').feature('sctr1').selection.set(sctrBnd);
model.component('comp1').physics('emw').feature('sctr1').label('Low-reflecting boundary, air top');

if ~strcmp(P.celltype,'rib')
    model.component('comp1').physics('emw').create('symp1', 'SymmetryPlane', 2);
    model.component('comp1').physics('emw').feature('symp1').selection.set(P.zEnd);
end
model.component('comp1').physics('emw').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('emw').feature('pc1').label('Periodic Condition x direction');
model.component('comp1').physics('emw').feature('pc2').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc2').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('emw').feature('pc2').label('Periodic Condition y direction');

% model.component('comp1').physics('emw').prop('components').set('components', 'inplane');

%% Add the solver and solver sequences 
study = model.study.create('std1');
std_param = study.create('param', 'Parametric');
std_eigv = study.create('eig', 'Eigenfrequency');
std_param.set('pname', 'k');
std_param.set('plistarr', kliststr);
std_param.set('punit', []);
std_eigv.set('neigsactive',true).set('neigs',nbands);
std_eigv.set('eigunit', 'THz');
std_eigv.set('shiftactive',true).set('shift',num2str(freq));

solv = model.sol.create('solv');
solv.study('std1');           % connect solver sequence to study node
solv.attach('std1');         % comes from .m saved from GUI - needed?

solv_stdstep = solv.create('st1', 'StudyStep'); % define study step, vars, solver node
solv_vars = solv.create('v1', 'Variables');
solv_eigv = solv.create('e1', 'Eigenvalue');
solv_eigv.create('d1','Direct');
study.label('Compile Equations: Eigenfrequency');
solv_vars.label('Dependent Variables 1.1');
solv_eigv.label('Eigenvalue Solver 1.1');

solv_eigv.set('transform','eigenfrequency');
% Units are given explicitly because the solver node does NOT inherit the
% study node's eigunit ('THz', set on std_eigv above): a bare num2str(freq)
% is read as Hz here, i.e. 1e12 times too small.
solv_eigv.set('shift', [num2str(freq),'[THz]']);
% eigref is the eigenvalue LINEARIZATION point - the frequency at which
% frequency-dependent quantities are evaluated - as distinct from shift, which
% only says where to search. It defaults to 0, and at 0 the Scattering
% condition's kinp (sctr1 below the physics section) evaluates
% sigma/(epsr*iomega*eps0) and divides by zero before the first iteration.
% The Sellmeier refractive index is frequency-dependent too and would likewise
% be evaluated at DC, far outside the fit's range of validity.
solv_eigv.set('eigref', [num2str(freq),'[THz]']);
solv_eigv.feature('dDef').label('Direct 2');
solv_eigv.feature('aDef').label('Advanced 1');
solv_eigv.feature('aDef').set('complexfun', true);
solv_eigv.feature('d1').label('Suggested Direct Solver (emw)');

% solv_eigv.set('transform', 'eigenfrequency');
% solv_eigv.set('control','std_eigv');
% solv_eigv.feature('dDef').set('linsolver', 'spooles');
% solv_eigv.feature('aDef').set('complexfun', 'off');

psolv = model.sol.create('psolv');
psolv.study('std1');
psolv.label('Parametric Solutions 2');

% add batch job configuration for parameter sweep
pbatch = model.batch.create('p1', 'Parametric');
pbatch_solseq = pbatch.create('so1', 'Solutionseq');
pbatch.study('std1');
pbatch.attach('std1');
pbatch.set('control', 'param');
pbatch.set('pname', 'k');
pbatch.set('plistarr', kliststr);
pbatch.set('punit', []);
pbatch.set('err', true);
pbatch_solseq.set('seq', 'solv');
pbatch_solseq.set('psol', 'psolv');
pbatch_solseq.set('param', kparamstr);

disp('Study, solver, sweep and batch nodes added');

%% Add Mesh
mesh_quality = meshSize;
mesh = model.mesh.create('mesh', 'geom1');
display(['Meshing with quality: ' num2str(mesh_quality)]);
mesh.autoMeshSize(mesh_quality).run;

%% Solve for bands
solv.runAll;
pbatch.run;

%% set up results node 
model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');
model.result.table.create('tbl7', 'Table');
model.result.table.create('tbl8', 'Table');
model.result.table.create('tbl9', 'Table');
model.result.table.create('tbl10', 'Table');
model.result.table.create('tbl11', 'Table');
model.result.table.create('tbl12', 'Table');
model.result.table.create('tbl13', 'Table');
model.result.table.create('tbl14', 'Table');


model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('gev2', 'EvalGlobal');
model.result.numerical('gev1').set('probetag', 'none');
model.result.numerical('gev2').set('data', 'dset2');
model.result.numerical('gev2').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup3D');
% model.result.create('pg2', 'PlotGroup2D');
% model.result.create('pg3', 'PlotGroup1D');
% model.result('pg1').create('surf1', 'Surface');
% model.result('pg2').set('data', 'dset2');
% model.result('pg2').create('surf1', 'Surface');
% model.result('pg3').set('data', 'dset2');
% model.result('pg3').create('glob1', 'Global');

model.result.numerical('gev1').label('Eigenfrequencies (emw)');
model.result.numerical('gev1').set('table', 'tbl13');
model.result.numerical('gev1').set('expr', {'emw.freq' 'emw.Qfactor'});
model.result.numerical('gev1').set('unit', {'THz' '1'});
model.result.numerical('gev1').set('descr', {'Frequency' 'Quality factor'});
model.result.numerical('gev2').label('Eigenfrequencies (emw) 1');
model.result.numerical('gev2').set('table', 'tbl14');
model.result.numerical('gev2').set('expr', {'emw.freq' 'emw.Qfactor'});
model.result.numerical('gev2').set('unit', {'THz' '1'});
model.result.numerical('gev2').set('descr', {'Frequency' 'Quality factor'});
model.result.numerical('gev1').setResult;
model.result.numerical('gev2').setResult;
% model.result('pg1').label('Electric Field (emw)');
% model.result('pg1').set('frametype', 'spatial');
% model.result('pg1').set('showlegendsmaxmin', true);
model.result('pg1').create('slc1', 'Slice');
model.result('pg1').label('Electric Field (emw)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('showlegendsmaxmin', true);
model.result('pg1').feature('slc1').set('quickplane', 'xy');
model.result('pg1').feature('slc1').set('quickzmethod', 'coord');
model.result('pg1').feature('slc1').set('quickz', 0);
model.result('pg1').feature('slc1').set('resolution', 'normal');

%% Save data and plots
if ~exist([P.datLoc,P.fileBase],'dir') && P.saveplots
    mkdir([P.datLoc,P.fileBase])
end

% extract solution info from parameter sweep
sols = mphsolutioninfo(model);
lambda_inds = find(strcmp(sols.psolv.mapheaders,'lambda'));
k_inds = find(strcmp(sols.psolv.mapheaders,'kx'));
inner_inds = find(strcmp(sols.psolv.mapheaders,'Inner'));
outer_inds = find(strcmp(sols.psolv.mapheaders,'Outer'));

% assemble solutions
for ki = 0:3*kpts-1
%     fem = mbfem;
    % assemble eigenvalues and eigenfrequencies
    lambda_ki = find(sols.psolv.map(:,outer_inds)==ki+1);
    fem.sol.lambda = sols.psolv.map(lambda_ki,lambda_inds);
    fem.sol.freqs = abs(fem.sol.lambda)/(2*pi);
    for nb = 1:nbands
        ds.F(ki+1,nb) = fem.sol.freqs(nb);
    end
end

% postprocess F and k
% append results from Gamma-point simulations to end of array
ds.F(end+1,1:nbands) = ds.F(1,1:nbands);

if P.bandStruct_2D
    % for 2D band structure
    ds.k_norm(end+1) = 3;
    ds.kx_norm(end+1) = ds.kx_norm(1);
    ds.ky_norm(end+1) = ds.ky_norm(1); 
else
    % for 1D band structures
    ds.kx_norm(end+1) = ds.kx_norm(1);
    % ky_norm is wrapped too, so it stays the same length as kx_norm on BOTH
    % branches. solveOpticalBands now builds the light line from
    % hypot(kx_norm,ky_norm) and would otherwise fail here on a size mismatch -
    % boomerang_optimize_sweep_diamond.m, rib_optimize_sweep.m and
    % rib_optimize_sweep_diamond.m all reach this branch with bandStruct_2D = 0.
    ds.ky_norm(end+1) = ds.ky_norm(1);
    ds.k_norm = ds.kx_norm;     
end

%% saving the mph files for debugging purposes 
if P.saveMPH
    path_mph = [P.datLoc,P.fileBase,'\bands_',txt_zsym,'.mph'];
    mphsave(model, path_mph);
end

ds.P = P;

end

% -------------------------------------------------------------------------

function bnd = airTopBoundary(model,P)
%AIRTOPBOUNDARY Boundaries of the air region's top face, normal parallel to z.
%
% Returns the indices to put the low-reflecting Scattering condition on. The
% boundary is found geometrically with bndindex, which takes a point in the plane
% plus a normal vector and returns every boundary coplanar with it - so
% [0 0 zTop] with normal [0 0 1] means "the face at height zTop whose normal is
% parallel to z". That is the same call buildSnowflakeUnitCell_optical.m:172
% already uses to set P.zEnd2, and it is why this survives a change in domain or
% boundary numbering that a literal index would not.
%
% P.zEnd2 wins when the builder already recorded it, so builders that locate the
% face themselves stay authoritative and are not second-guessed here.

if isfield(P,'zEnd2') && ~isempty(P.zEnd2)
    bnd = P.zEnd2;
    return
end

zTop = airTopZ(P);
geom = model.component('comp1').geom('geom1');
bnd  = bndindex(geom,[0 0 zTop],[0 0 1]);

if isempty(bnd)
    error('runOpticalBand_3D:noAirTopBoundary', ...
        ['No boundary with a z-parallel normal was found at z = %g nm, so the ' ...
         'low-reflecting condition has nothing to attach to. Either the air ' ...
         'region was not built (check P.add_airDisk in the geometry builder) or ' ...
         'its top does not sit where P implies.'], zTop*1e9);
end

% A Scattering condition belongs on an EXTERIOR boundary. If the builder also
% stacks an absorbing slab above the air - buildSnowflakeUnitCell_optical puts
% its PML work plane at exactly this height - then the face found here is
% interior, and the condition should go on the outermost face instead. Flagged
% rather than auto-corrected, because which face is outermost is a property of
% the builder, not of P.
fprintf('  low-reflecting boundary at z = %.1f nm, %d face(s)\n', ...
    zTop*1e9, numel(bnd));
end

% -------------------------------------------------------------------------

function zTop = airTopZ(P)
%AIRTOPZ z_max of the air disk [m], matching buildBoomerangUnitCell's convention.
%
% The air disk is extruded from the base of the remaining solid upward by its
% height, so its top is base + height. The base tracks the z-symmetry state: the
% slab spans -th/2 to +th/2, but a nonzero P.mbevenz means the builder subtracted
% everything below z = 0, so the air starts at 0 in that case and at -th/2
% otherwise. Keep this in step with the air disk block in
% buildBoomerangUnitCell.m if that convention ever changes.

if ~(isfield(P,'add_airDisk') && P.add_airDisk)
    error('runOpticalBand_3D:airDiskDisabled', ...
        ['The low-reflecting boundary is placed on the top of the air disk, but ' ...
         'P.add_airDisk is off or absent so no air disk exists. Enable it in the ' ...
         'test script, or set P.zEnd2 in the geometry builder to name the face ' ...
         'directly.']);
end

if isfield(P,'airDiskHeight')
    airDiskHeight = P.airDiskHeight;
elseif isfield(P,'airDiskH')      % alias used by the rest of this directory
    airDiskHeight = P.airDiskH;
else
    error('runOpticalBand_3D:noAirDiskHeight', ...
        'P.add_airDisk is on but neither P.airDiskHeight nor P.airDiskH is set.');
end
validateattributes(airDiskHeight,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'air disk height');

if isfield(P,'mbevenz') && abs(P.mbevenz)
    zAirBase = 0;
else
    zAirBase = -P.th/2;
end
zTop = zAirBase + airDiskHeight;
end
