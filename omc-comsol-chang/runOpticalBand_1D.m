function ds = runOpticalBand_1D(P)

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

ModelUtil.clear();
clear model     
model = ModelUtil.create('model');

beamMat = P.beamMat;
a = P.a;

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

%% Define k-points for sweep over wavevectors (1D band structure)
% adapted from phononic crystal model on COMSOL

% parameter node for COMSOL model
% k runs from 0 to 3: 0-->1 for Gamma-X
model.param.set('k', '1');
model.param.set('a', [num2str(a),'[m]']);

if strcmp(P.unitcell,'hexagonal')
    model.param.set('kx', '(pi/a)*k*(sqrt(3)/2)');
    model.param.set('ky', '(pi/a)*k*(-1/2)');
    
    for ki = 0:kpts
        ds.k_norm(ki+1,1) = ki/kpts;
        ds.kx_norm(ki+1,1) = ((sqrt(3)/2)*ki/kpts)*(ki<kpts);                  % Gamma-X
        ds.ky_norm(ki+1,1) = 0;                          % Gamma-X
                            
    end
    
    % compile expressions for input to COMSOL model
    kliststr = ['range(0,1/',num2str(kpts),',1)'];
    % kliststr = ['range(0,1/',num2str(kpts),',1)'];
    for ki = 1:kpts
        kparamstr{ki} = ['"k", "',num2str((ki-1)/kpts),'"'];
    end
else
    model.param.set('kx', 'pi/a*k');
    model.param.set('ky', '0');
    
    for ki = 0:kpts-1
        ds.k_norm(ki+1,1) = ki/kpts;
        ds.kx_norm(ki+1,1) = (ki/kpts)*(ki<kpts);                  % Gamma-X
        ds.ky_norm(ki+1,1) = 0;                         % Gamma-X
    end
    % compile expressions for input to COMSOL model
    kliststr = ['range(0,1/',num2str(kpts),',1)'];
    for ki = 1:kpts
        kparamstr{ki} = ['"k", "',num2str((ki-1)/kpts),'"'];
    end
end 

%% Set up the geometry
if strcmp(P.celltype,'cross')
    [model,P] = DrawCrossUnitCell(model,P);
elseif strcmp(P.celltype,'boomerang')
    [model,P] = buildBoomerangUnitCellStrip_v2(model,P);
elseif strcmp(P.celltype,'boomerang_strip_v2')
    [model,P] = buildBoomerangUnitCellStrip_v2(model,P);
elseif strcmp(P.celltype,'boomerang_strip')
    [model,P] = buildBoomerangUnitCellStrip(model,P);
elseif strcmp(P.celltype,'hole_strip')
    [model,P] = buildHoleStrip_3D(model,P);
elseif strcmp(P.celltype,'hole_strip_wvg')
    [model,P] = buildHoleStrip_withWg_2D(model,P);
elseif strcmp(P.celltype,'Snowflake_strip_2d')
    [model,P] = buildSnowflakeStrip_2D(model,P);
elseif strcmp(P.celltype,'rib')
    [model,P] = buildRibUnitCell_LN(model,P);
elseif strcmp(P.celltype,'hole')
    [model,P] = buildHoleUnitCell(model,P);
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
model.component('comp1').geom('geom1').run;
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').selection.all;
model.component('comp1').material.create('mat1', 'Common');
if strcmp(P.celltype,'boomerang_strip_v2')
    model.component('comp1').material('mat1').selection.named([P.ucellname,'_beamSel']);
else
    model.component('comp1').material('mat1').selection.set([1]);
end
if strcmp(P.celltype,'rib')
    model.component('comp1').material('mat1').label('Air');
    model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'10e-12' '0' '0' '0' '10e-12' '0' '0' '0' '10e-12'});
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    if strcmp(P.beamMat,'LN')
    model.component('comp1').material('mat2').label('LiNbO3 (Lithium niobate)');
    model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', '');
    model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', '');
    model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func.create('an1', 'Analytic');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func.create('an2', 'Analytic');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').label('Sellmeyer formula - extraordinary refractive index');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').set('funcname', 'neref_sel');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').set('expr', '(1 + ((lambda0^(2))*2.9804e12)/((lambda0^(2))*1e12 - 0.02047) + ((lambda0^(2))*0.5981e12)/((lambda0^(2))*1e12 - 0.0666) + ((lambda0^(2))*8.9543e12)/((lambda0^(2))*1e12 - 416.08))^(1./2.)');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').set('args', {'lambda0'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').set('argunit', {'m'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an1').set('plotargs', {'lambda0' '0.4[um]' '5[um]'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').label('Sellmeyer formula - ordinary refractive index');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').set('funcname', 'noref_sel');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').set('expr', '(1 + ((lambda0^(2))*2.6734e12)/((lambda0^(2))*1e12 - 0.01764) + ((lambda0^(2))*1.22901e12)/((lambda0^(2))*1e12 - 0.05914) + ((lambda0^(2))*12.614e12)/((lambda0^(2))*1e12 - 474.60))^(1./2.)');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').set('args', {'lambda0'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').set('argunit', {'m'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').func('an2').set('plotargs', {'lambda0' '0.4[um]' '5[um]'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'neref_sel(lbd0)' '0' '0' '0' 'noref_sel(lbd0)' '0' '0' '0' 'noref_sel(lbd0)'});
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').addInput('frequency');
    elseif strcmp(P.beamMat,'diamond')
        model.component('comp1').material('mat2').label('diamond');
    model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
        model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'2.4028' '0' '0' '0' '2.4022' '0' '0' '0' '2.4028'});
        model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0.001' '0' '0' '0' '0.001' '0' '0' '0' '0.001'});
        model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'5.7734' '0' '0' '0' '5.7734' '0' '0' '0' '5.7734'});
    end
else
    model.component('comp1').material('mat2').label('Air');
    model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'10e-12' '0' '0' '0' '10e-12' '0' '0' '0' '10e-12'});
    model.component('comp1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
    if strcmp(P.beamMat,'diamond')
        model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.4028' '0' '0' '0' '2.4022' '0' '0' '0' '2.4028'});
        model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0.001' '0' '0' '0' '0.001' '0' '0' '0' '0.001'});
        model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'5.7734' '0' '0' '0' '5.7734' '0' '0' '0' '5.7734'});
    elseif strcmp(P.beamMat,'SiC')
        model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.5' '0' '0' '0' '2.5' '0' '0' '0' '2.5'});
    else 
        model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'3.5' '0' '0' '0' '3.5' '0' '0' '0' '3.5'});
    end
end

%% setup the physics and the boundary conditions
% model.component('comp1').coordSystem.create('pml1', 'PML');
% model.component('comp1').coordSystem('pml1').selection.set([3]);
model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('emw').feature('pc1').selection.named('geom1_xboundaries_bnd');
if P.bandStructureDim == 1
    % for TE mode (PEC in the y direction)
    model.component('comp1').physics('emw').create('sympy', 'SymmetryPlane', 2);
    model.component('comp1').physics('emw').feature('sympy').selection.named('geom1_yboundaries_bnd');
    model.component('comp1').physics('emw').feature('sympy').set('Symmetry_type', 'pec');
else
    model.component('comp1').physics('emw').create('pc2', 'PeriodicCondition', 2);
    model.component('comp1').physics('emw').feature('pc2').selection.named('geom1_yboundaries_bnd');
end
if ~strcmp(P.celltype,'rib') && P.mbevenz   
    model.component('comp1').physics('emw').create('symp1', 'SymmetryPlane', 2);
    model.component('comp1').physics('emw').feature('symp1').selection.named('geom1_ZsymSel');
end
model.component('comp1').physics('emw').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('emw').feature('pc1').label('Periodic Condition x direction');
if P.bandStructureDim ~= 1
    model.component('comp1').physics('emw').feature('pc2').set('PeriodicType', 'Floquet');
    model.component('comp1').physics('emw').feature('pc2').set('kFloquet', {'kx'; 'ky'; '0'});
    model.component('comp1').physics('emw').feature('pc2').label('Periodic Condition y direction');
end
% add the scattering boundary condition 
model.component('comp1').physics('emw').create('sctr1', 'Scattering', 2);
if strcmp(P.celltype,'boomerang_strip_v2')
    model.component('comp1').physics('emw').feature('sctr1').selection.set([7]);
else
    model.component('comp1').physics('emw').feature('sctr1').selection.set([3 12]);
end

%% Add Mesh
mesh_quality = meshSize;
mesh = model.mesh.create('mesh', 'geom1');
display(['Meshing with quality: ' num2str(mesh_quality)]);
mesh.autoMeshSize(mesh_quality).run;

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

% sol1 = model.sol.create('sol1');
% sol1.study('std1');
% sol1.attach('std1');
% sol1.create('st1', 'StudyStep');
% sol1.create('v1', 'Variables');
% sol1.create('e1', 'Eigenvalue');
% sol1.feature('e1').create('d1', 'Direct');
% sol1.feature('e1').create('i1', 'Iterative');
% sol1.feature('e1').set('shift', num2str(freq));
% sol1.feature('e1').feature('i1').create('mg1', 'Multigrid');
% sol1.feature('e1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
% sol1.feature('e1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');
% sol1.feature('e1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');

solv_stdstep = solv.create('st1', 'StudyStep'); % define study step, vars, solver node
solv_vars = solv.create('v1', 'Variables');
solv_eigv = solv.create('e1', 'Eigenvalue');
% solv_eigv.set('neigsactive',true).set('neigs',nbands);
model.sol('solv').feature('e1').set('shift', [num2str(freq)]);
model.sol('solv').feature('e1').set('neigs', nbands);
solv_eigv.create('d1','Direct');
study.label('Compile Equations: Eigenfrequency');
solv_vars.label('Dependent Variables 1.1');
solv_eigv.label('Eigenvalue Solver 1.1');

solv_eigv.set('transform','eigenfrequency');
solv_eigv.set('shift',num2str(freq));
solv_eigv.feature('dDef').label('Direct 2');
solv_eigv.feature('aDef').label('Advanced 1');
solv_eigv.feature('aDef').set('complexfun', true);
solv_eigv.feature('d1').label('Suggested Direct Solver (emw)');

% solv_eigv.set('transform', 'eigenfrequency');
% solv_eigv.set('control','std_eigv');
% solv_eigv.feature('dDef').set('linsolver', 'spooles');
% solv_eigv.feature('aDef').set('complexfun', 'off');

% psolv = model.sol.create('psolv');
% psolv.study('std1');
% psolv.label('Parametric Solutions 2');

% add batch job configuration for parameter sweep
pbatch = model.batch.create('psolv', 'Parametric');
pbatch_solseq = pbatch.create('so1', 'Solutionseq');
pbatch.study('std1');
pbatch.attach('std1');
pbatch.set('control', 'param');
pbatch.set('pname', 'k');
pbatch.set('plistarr', kliststr);
pbatch.set('punit', []);
pbatch.set('err', true);
pbatch_solseq.set('seq', 'solv');
% pbatch_solseq.set('psol', 'psolv');
pbatch_solseq.set('param', kparamstr);
model.study('std1').feature('param').set('pname', {'k'});
model.study('std1').feature('param').set('plistarr', {kliststr});
model.study('std1').feature('param').set('punit', {''});
model.study('std1').feature('eig').set('neigs', nbands);
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('eigunit', 'THz');
model.study('std1').feature('eig').set('shift', num2str(P.optical_freq));

model.sol('solv').feature('e1').set('control', 'eig');
model.sol('solv').feature('e1').set('transeigref', true);
model.sol('solv').feature('e1').set('eigref', [num2str(freq)]);
% model.sol('solv').feature('e1').set('eigunit', 'THz');

disp('Study, solver, sweep and batch nodes added');
% study = model.study.create('std1');
% std_param = study.create('param', 'Parametric');
% std_eigv = study.create('eig', 'Eigenfrequency');
% std_param.set('pname', 'k');
% std_param.set('plistarr', kliststr);
% std_param.set('punit', []);
% 
% solv = model.sol.create('solv');
% solv.study('std1');           % connect solver sequence to study node
% solv.attach('std1');         % comes from .m saved from GUI - needed?
% 
% model.sol.create('sol1');
% model.sol('sol1').study('std1');
% model.sol('sol1').attach('std1');
% model.sol('sol1').create('st1', 'StudyStep');
% model.sol('sol1').create('v1', 'Variables');
% model.sol('sol1').create('e1', 'Eigenvalue');
% model.sol('sol1').feature('e1').create('d1', 'Direct');
% model.sol('sol1').feature('e1').create('i1', 'Iterative');
% model.sol('sol1').feature('e1').feature('i1').create('mg1', 'Multigrid');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
% model.sol.create('sol2');
% model.sol('sol2').study('std1');
% model.sol('sol2').label('Parametric Solutions 1');
% 
% 
% psolv = model.sol.create('psolv');
% psolv.study('std1');
% psolv.label('Parametric Solutions 2');
% 
% % add batch job configuration for parameter sweep
% pbatch = model.batch.create('p1', 'Parametric');
% pbatch_solseq = model.batch('p1').create('so1', 'Solutionseq');
% model.batch('p1').study('std1');
% pbatch.attach('std1');
% pbatch.set('control', 'param');
% pbatch.set('pname', 'k');
% pbatch.set('plistarr', kliststr);
% pbatch.set('punit', '');
% pbatch.set('err', true);
% pbatch_solseq.set('seq', 'solv');
% pbatch_solseq.set('psol', 'psolv');
% pbatch_solseq.set('param', kparamstr);
% 
% model.study('std1').feature('param').set('pname', {'k'});
% model.study('std1').feature('param').set('plistarr', {kliststr});
% model.study('std1').feature('param').set('punit', {''});
% model.study('std1').feature('eig').set('neigs', nbands);
% model.study('std1').feature('eig').set('neigsactive', true);
% model.study('std1').feature('eig').set('eigunit', 'THz');
% model.study('std1').feature('eig').set('shift', num2str(P.optical_freq));
% 
% model.sol('sol1').attach('std1');
% model.sol('sol1').feature('st1').label('Compile Equations: Eigenfrequency');
% model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
% model.sol('sol1').feature('e1').label('Eigenvalue Solver 1.1');
% model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
% model.sol('sol1').feature('e1').set('neigs', nbands);
% model.sol('sol1').feature('e1').set('shift',  num2str(P.optical_freq));
% model.sol('sol1').feature('e1').set('eigref', num2str(P.optical_freq));
% model.sol('sol1').feature('e1').feature('dDef').label('Direct 2');
% model.sol('sol1').feature('e1').feature('aDef').label('Advanced 1');
% model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
% model.sol('sol1').feature('e1').feature('d1').active(true);
% model.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (emw)');
% model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'pardiso');
% model.sol('sol1').feature('e1').feature('i1').label('Suggested Iterative Solver (emw)');
% model.sol('sol1').feature('e1').feature('i1').set('itrestart', 300);
% model.sol('sol1').feature('e1').feature('i1').set('prefuntype', 'right');
% model.sol('sol1').feature('e1').feature('i1').feature('ilDef').label('Incomplete LU 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').label('Multigrid 1.1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').set('iter', 1);
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').feature('sv1').label('SOR Vector 1.1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').feature('sv1').label('SOR Vector 1.1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
% model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
% % 
% % 
% % std_param.set('pname', 'k');
% % std_param.set('plistarr', kliststr);
% % std_param.set('punit', '');
% % model.study('std1').feature('eig').set('eigunit', 'THz');
% % std_eigv.set('neigsactive',true).set('neigs',nbands);
% % model.study('std1').feature('eig').set('shift',num2str(freq));
% % 
% % solv = model.sol.create('solv');
% % solv.study('std1');           % connect solver sequence to study node
% % solv.attach('std1');         % comes from .m saved from GUI - needed?
% % 
% % solv_stdstep = solv.create('st1', 'StudyStep'); % define study step, vars, solver node
% % solv_vars = solv.create('v1', 'Variables');
% % solv_eigv = solv.create('e1', 'Eigenvalue');
% % solv_eigv.create('d1','Direct');
% % study.label('Compile Equations: Eigenfrequency');
% % solv_vars.label('Dependent Variables 1.1');
% % solv_eigv.label('Eigenvalue Solver 1.1');
% % 
% % solv_eigv.set('transform','eigenfrequency');
% % solv_eigv.set('shift',num2str(freq));
% % solv_eigv.feature('dDef').label('Direct 2');
% % solv_eigv.feature('aDef').label('Advanced 1');
% % solv_eigv.feature('aDef').set('complexfun', true);
% % solv_eigv.feature('d1').label('Suggested Direct Solver (ewfd)');
% % 
% % 
% % % solv_eigv.set('transform', 'eigenfrequency');
% % % solv_eigv.set('control','std_eigv');
% % % solv_eigv.feature('dDef').set('linsolver', 'spooles');
% % % solv_eigv.feature('aDef').set('complexfun', 'off');
% % 
% % psolv = model.sol.create('psolv');
% % psolv.study('std1');
% % psolv.label('Parametric Solutions 2');
% % 
% % add batch job configuration for parameter sweep
% 
% disp('Study, solver, sweep and batch nodes added');


%% Solve for bands
solv.runAll;
% model.sol('sol1').run;
pbatch.run;

%% set up results node 
% model.result.table.create('tbl1', 'Table');
% model.result.table.create('tbl2', 'Table');
% model.result.table.create('tbl3', 'Table');
% model.result.table.create('tbl4', 'Table');
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
% visualizing the simulation results 
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').set('data', 'dset2');
model.result('pg1').create('mslc1', 'Multislice');
model.result('pg1').feature('mslc1').create('filt1', 'Filter');
model.result('pg1').feature('mslc1').feature('filt1').set('expr', '!isScalingSystemDomain');
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
lambda_inds = find(strcmp(sols.sol1.mapheaders,'lambda'));
k_inds = find(strcmp(sols.solv.mapheaders,'k'));
inner_inds = find(strcmp(sols.sol1.mapheaders,'Inner'));
outer_inds = find(strcmp(sols.sol1.mapheaders,'Outer'));

% assemble solutions
for ki = 1:kpts
%     fem = mbfem;
    % assemble eigenvalues and eigenfrequencies
    lambda_ki = find(sols.sol1.map(:,outer_inds)==ki+1);
    fem.sol.lambda = sols.sol1.map(lambda_ki,lambda_inds);
    fem.sol.freqs = abs(fem.sol.lambda)/(2*pi);
    for nb = 1:nbands
        ds.F(ki+1,nb) = fem.sol.freqs(nb);
    end
end

% postprocess F and k
% append results from Gamma-point simulations to end of array
% ds.F(end+1,1:nbands) = ds.F(1,1:nbands);

% ds.k_norm(end+1) = ds.k_norm(1);
ds.kx_norm(end+1) = 1;
ds.k_norm = ds.kx_norm;     % for 1D band structures

%% saving the mph files for debugging purposes 
if P.saveMPH
    path_mph = [P.datLoc,P.fileBase,'\bands_',txt_zsym,'.mph'];
    mphsave(model, path_mph);
end

%% export the raw data 
if P.saveRawData
    model.result.export.create('tbl_exp', 'Table');
    model.result.export('tbl_exp').label('bandStruct');
    model.result.export('tbl_exp').set('table', 'tbl14');
    model.result.export('tbl_exp').set('filename', '.\bandStruct_data\BandStruct.txt');
    model.result.export('tbl_exp').set('header', false);
    model.result.export('tbl_exp').set('notation', 'scientific');
    model.result.export('tbl_exp').run;
end 

ds.P = P; 

end 
