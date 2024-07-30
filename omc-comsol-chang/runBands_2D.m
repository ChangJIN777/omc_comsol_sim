function ds = runBands_2D(P)
%RUNBANDS Summary of this function goes here
%   Detailed explanation goes here
% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

ModelUtil.clear();
clear model     

% create COMSOL model named 'model' from which COMSOL methods can be called, 
% e.g. model.save
model = ModelUtil.create('model');
% model.component.create('comp1', true);
% extract parameters from P
beamMat = P.beamMat;
a = P.a;

kpts = P.kpts;
nbands = P.nbands;
freq = P.freq;
evenz = P.mbevenz;
eveny = P.mbeveny;
max_dof = P.max_dof;
meshSize = P.meshSize;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric


% prefixes for filenames
if evenz == 1
    txt_zsym = 'evenz';
    txt_sol = 'even z';
elseif evenz == -1
    txt_zsym = 'oddz';
    txt_sol = 'odd z';
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
        ds.ky_norm(ki+1,1) = ((-1/2)*(ki<kpts)+...                          % Gamma-X
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
    [model,P] = buildLowerBoomerangUnitCell(model,P);
elseif strcmp(P.celltype,'snowflake')
    [model,P] = buildSnowflakeUnitCell(model,P);
end

if P.plotgeom
    figure;
    mphgeom(model)
    pathFig = [P.datLoc,'\',P.fileBase,'_geom'];
    saveas(gcf,[pathFig,'.fig']);
    saveas(gcf,[pathFig,'.png']);
end

%% Define material and properties
matTags = mphmodel(model.material);
P = LoadMaterialParams(P);
% create material nodes
if ~isfield(matTags,beamMat)
    bMat = model.material.create(beamMat);
else
    bMat = model.material(beamMat);
end
bMat.label(beamMat);
bMat_def = bMat.propertyGroup('def');

% density
mbfem.rho = P.rho;
bMat_def.set('density', mbfem.rho);

% rotate crystal from [100] direction, about normal to {100}
if (isfield(P,'rxtal') && isfield(P,'D'))
    [~,D] = RotateXtalTensor(P.D,P.rxtal);
    P.D = D;
    ds.P = P;      
end

% define anisotropic material (and generate full 6x6 elasticity matrix), 
% or define material based on Young's modulus and Poisson's ratio
if isfield(P,'D')
    if length(P.D) < 21
        error('Insufficient number of stiffness constants specified');
    end
    bMat_aniso = bMat.propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
    bMat_aniso.set('DVo', P.D);
    
elseif (isfield(P,'E') && isfield(P,'nu')) 
    mbfem.nu = P.nu;
    mbfem.E = P.E;
    bMat_def.set('youngsmodulus', mbfem.E);
    bMat_def.set('poissonsratio', mbfem.nu);
else
    error('No stiffness type specified');
end

% jagged_array = NET.createArray('System.Int32[]',1);
% domain_num = [1];
% jagged_array(1) = domain_num;
mbfem.b_domind = 1;
bMat.selection.geom('geom1', 3);
model.component('comp1').geom('geom1').run;
bMat.selection.set(mbfem.b_domind);


%% Setup the physics and boundary conditions
smech = model.physics.create('smech', 'SolidMechanics', 'geom1');
if (isfield(P,'rxtal') && isfield(P,'D'))
    smech.feature('lemm1').set('SolidModel', 'Anisotropic');
    smech.feature('lemm1').set('AnisotropicOption', 'AnisotropicVo');
    smech.feature('lemm1').set('DVo_mat', 'from_mat');
elseif (isfield(P,'E') && isfield(P,'nu')) 
    smech.feature('lemm1').set('SolidModel', 'Isotropic');
    smech.feature('lemm1').set('IsotropicOption', 'Enu');
    smech.feature('lemm1').set('E_mat', 'from_mat');
    smech.feature('lemm1').set('nu_mat', 'from_mat');
end
smech.feature('lemm1').set('rho_mat', 'from_mat');
smech.selection.set(mbfem.b_domind);
% lem = linear elastic material

% boundary conditions
% all BCs set to free by default
clear bnds
bnds.pbc_inds = [];
if P.fixed_bc
    bnds.fixed_inds = [];
    % fixed BCs for xz planes at y = +/- w/2
    fixedBCs = smech.create('fix1', 'Fixed', 2);
    fixedY1 = P.yEnd1;
    fixedY2 = P.yEnd2;
    fixedYinds = [fixedY1 fixedY2];
    bnds.fixed_inds(end+1:end+length(fixedYinds)) = fixedYinds;
    
    % % fixed BCs for xz planes at x = +/- w/2
    % fixedX1 = P.xEnd1;
    % fixedX2 = P.xEnd2;
    % % fixedXinds = [fixedX1 fixedX2];
    % fixedXinds = [];
    % bnds.fixed_inds(end+1:end+length(fixedXinds)) = fixedXinds;
    
    % set the fixed BCs 
    fixedBCs.selection.set(fixedXinds);
end 

% periodic BCs for yz planes at x = +/- a/2
pbcX = smech.create('pbcX', 'PeriodicCondition', 2);
pbcX.label('Periodic BC, x-direction');
pbcX.set('PeriodicType', 'Floquet');
pbcX1 = P.xEnd1;
pbcX2 = P.xEnd2;
pbcXinds = [pbcX1 pbcX2];
bnds.pbc_inds(end+1:end+length(pbcXinds)) = pbcXinds;
pbcX.selection.set(pbcXinds);
pbcX.set('manualDestinationSelection', true);
pbcX.selection('destinationDomains').set(pbcX2);
pbcX.set('kFloquet', {'kx'; 'ky'; '0'});        %initialize Floquet vector (2D band structures)

% periodic BCs for xz planes at y = +/- w/2 (used to simulate 2D band structures)
pbcY = smech.create('pbcY', 'PeriodicCondition', 2);
pbcY.label('Periodic BC, y-direction');
pbcY.set('PeriodicType', 'Floquet');
pbcY1 = P.yEnd1;
pbcY2 = P.yEnd2;
pbcYinds = [pbcY1 pbcY2];
bnds.pbc_inds(end+1:end+length(pbcYinds)) = pbcYinds;
pbcY.selection.set(pbcYinds);
pbcY.set('manualDestinationSelection', true);
pbcY.selection('destinationDomains').set(pbcY2);
pbcY.set('kFloquet', {'kx'; 'ky'; '0'});        %initialize Floquet vector

% symmetric BC for xy plane containing point (0,0,0)
symBCs = smech.create('symBCs', 'SymmetrySolid', 2);
symBCs.label('Symmetric BC');
bnds.sym_inds = [];
if evenz == 1
    bndinds =P.zEnd;
    bnds.sym_inds(end+1:end+length(bndinds)) = bndinds;
    clear bndinds
end
if (~isempty(bnds.sym_inds))
    symBCs.selection.set(bnds.sym_inds);
end

% anti-symmetric BC for xy plane containing point (0,0,0)
asymBCs = smech.create('asymBCs', 'Antisymmetry', 2);
asymBCs.label('Anti-symmetric BC');
bnds.asym_inds = [];
if evenz == -1 
    bndinds = P.zEnd;
    bnds.asym_inds(end+1:end+length(bndinds)) = bndinds;
    clear bndinds
end
if (~isempty(bnds.asym_inds))
    asymBCs.selection.set(bnds.asym_inds);
end

mbfem.bnds = bnds;
disp('Solid Mechanics added - boundary conditions done');

%% Add the solver and the solver sequences
% add parametric and eigenfrequency study
study = model.study.create('study');
std_param = study.create('std_param', 'Parametric');
std_param.set('pname', 'k');
std_param.set('plistarr', kliststr);
std_param.set('punit', '');
std_eigv = study.create('std_eigv','Eigenfrequency');
std_eigv.set('neigsactive',true).set('neigs',nbands);
std_eigv.set('shiftactive',true).set('shift',num2str(freq));

solv = model.sol.create('solv');
solv.study('study');           % connect solver sequence to study node
solv.attach('study');         % comes from .m saved from GUI - needed?
solv_stdstep = solv.create('solv_stdstep', 'StudyStep'); % define study step, vars, solver node
solv_stdstep.set('study','study').set('studystep','std_eigv');
solv_vars = solv.create('solv_vars', 'Variables');
solv_vars.set('control','std_eigv');
solv_eigv = solv.create('solv_eigv', 'Eigenvalue');
solv_eigv.set('transform', 'eigenfrequency');
solv_eigv.set('control','std_eigv');
solv_eigv.set('eigref',num2str(freq));
solv_eigv.feature('dDef').set('linsolver', 'spooles');
solv_eigv.feature('aDef').set('complexfun', 'off');

psolv = model.sol.create('psolv');
psolv.study('study');

% add batch job configuration for parameter sweep
pbatch = model.batch.create('pbatch', 'Parametric');
pbatch_solseq = pbatch.create('pbatch_solseq', 'Solutionseq');
pbatch.study('study');
pbatch.attach('study');
pbatch.set('pname', 'k');
pbatch.set('plistarr', kliststr);
pbatch.set('punit', '');
pbatch.set('err', true);
pbatch.set('control', 'std_param');
pbatch_solseq.set('psol', 'psolv');
pbatch_solseq.set('param', kparamstr);
pbatch_solseq.set('seq', 'solv');

disp('Study, solver, sweep and batch nodes added');

%% Add Mesh
mesh_quality = meshSize;
mesh = model.mesh.create('mesh', 'geom1');
display(['Meshing with quality: ' num2str(mesh_quality)]);
mesh.autoMeshSize(mesh_quality).run;

%% Adjust mesh if max DOFs exceeded
mesh_ok = 0;
% Xmesh info - get degrees of freedom
mbfem.xmesh = mphxmeshinfo(model, 'soltag', 'solv', ...
                                   'studysteptag', 'solv_stdstep');
dofs = mbfem.xmesh.ndofs;
display(['Estimated no. of DoFs: ' num2str(dofs)]);

while (~mesh_ok) && (mesh_quality < 10)
    while (dofs >= max_dof) && (mesh_quality < 10)
        mesh_quality = mesh_quality + 1;
        display(['Meshing with quality: ' num2str(mesh_quality)]);
        mesh.autoMeshSize(mesh_quality).run;
        mbfem.xmesh = mphxmeshinfo(model, 'soltag', 'solv', ...
                                           'studysteptag', 'solv_stdstep');
        dofs = mbfem.xmesh.ndofs;
        display(['Estimated no. of DoFs: ' num2str(dofs)]);
    end
    mesh_ok = 1;
end
mphsave('test_geom')
mbfem.mbMesh = mesh_quality;

% debugging 
% mphlaunch(model);
%% Solve for bands
solv.runAll;
pbatch.run;

%% set up results node 

model.result.dataset('dset1').tag('dset');
dset = model.result.dataset('dset');
dset.set('solution', 'solv');
model.result.dataset('dset2').tag('pdset');
pdset = model.result.dataset('pdset');
pdset.set('solution', 'psolv');

if P.saveplots
    % % create sector dataset
    % xySymFac = 2^abs(eveny);
    % zSym = abs(evenz)*strcmp(P.xsect,'rect');
    % zSymFac = 2^zSym;
    % 
    % dsetTags = mphmodel(model.result.dataset);
    % 
    % if ~isfield(dsetTags,'pdset_sec')
    %     pdset_sec = model.result.dataset.create('pdset_sec', 'Sector3D');
    % else
    %     pdset_sec = model.result.dataset('pdset_sec');
    % end

    % pdset_sec.label('M Sol Full Beam XY');
    % pdset_sec.set('data', 'pdset');
    % pdset_sec.set('method', 'twopoint');
    % pdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
    % pdset_sec.setIndex('genpoints', '0', 0, 1);
    % pdset_sec.setIndex('genpoints', '0', 0, 2);
    % pdset_sec.setIndex('genpoints', '0', 1, 0);
    % pdset_sec.setIndex('genpoints', '0', 1, 1);
    % pdset_sec.setIndex('genpoints', '1', 1, 2);
    % if xySymFac
    %     pdset_sec.set('sectors', xySymFac);
    %     pdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
    %     pdset_sec.set('reflaxis', {'0' '0' '1'});   % reflection axis
    % end

    model.result.dataset.create('sec1', 'Sector3D');
    model.result.dataset('sec1').set('trans', 'rotrefl');
    model.result.dataset('sec1').set('pddir', {'1' '0' '0'});
    model.result.dataset('sec1').set('reflaxis', {'0' '1' '0'});
    model.result.dataset('sec1').set('data', 'pdset');
    if evenz == -1
        model.result.dataset('sec1').set('reflinv', 'on');
    end
    % if eveny==-1
    %     pdset_sec.set('rotinv', 'on');  % odd symmetry about x-plane
    %     pdset_sec.set('reflinv', 'on');
    % end
    % pdsetPlot = 'pdset_sec';

    % if zSym
    %     if ~isfield(dsetTags,'pdset_secZ')
    %         pdset_sec = model.result.dataset.create('pdset_secZ', 'Sector3D');
    %     else
    %         pdset_sec = model.result.dataset('pdset_secZ');
    %     end
    %     pdset_sec.label('M Sol Full Beam XYZ');
    %     pdset_sec.set('data', 'pdset_sec');
    %     pdset_sec.set('method', 'twopoint');
    %     pdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
    %     pdset_sec.setIndex('genpoints', '0', 0, 1);
    %     pdset_sec.setIndex('genpoints', '0', 0, 2);
    %     pdset_sec.setIndex('genpoints', '1', 1, 0);
    %     pdset_sec.setIndex('genpoints', '0', 1, 1);
    %     pdset_sec.setIndex('genpoints', '0', 1, 2);
    %     pdset_sec.set('sectors', zSymFac);
    %     pdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
    %     pdset_sec.set('reflaxis', {'0' '0' '1'});   % reflection axis
    %     if evenz==-1
    %         pdset_sec.set('rotinv', 'on');  % odd symmetry about z-plane
    %         pdset_sec.set('reflinv', 'on');
    %     end
    %     pdsetPlot = 'pdset_secZ';
    % 
    % end

    % if evenz == -1
    %     model.result.dataset('sec1').set('reflinv', 'on');
    % end
    
    % create 3D plot group for displacement field
    dispplot = model.result.create('dispplot', 'PlotGroup3D');
    dispplot.set('data','sec1');
    dispplot.set('solrepresentation', 'solutioninfo');
    dispplot.set('titletype', 'none');
    dispplot_vol = dispplot.create('dispplot_vol', 'Volume');
    dispplot_vol.set('rangedataactive','on').set('rangecoloractive','on');
    dispplot_def = dispplot_vol.create('dispplot_def', 'Deform');
    dispplot_vol.set('data', 'parent');
    dispplot.run;

    % % create 3D plot group for strain profile
    % strplot = model.result.create('strplot', 'PlotGroup3D');
    % strplot.set('data',pdsetPlot);
    % strplot.set('solrepresentation', 'solutioninfo');
    % strplot.set('titletype', 'none');
    % strplot_slc = strplot.create('strplot_slc', 'Slice');
    % strplot_slc.set('planetype', 'quick');
    % strplot_slc.set('quickxmethod', 'coord').set('quickx', (1-holeatedge/2)*a);
    % strplot_slc.set('rangedataactive','on').set('rangecoloractive','on');
    % strplot_slc.set('data', 'parent');
    % strplot.run;
    
end

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


    % if P.bandStruct_2D
    %     if strcmp(P.unitcell,'hexagonal')
    % 
    %         ds.k_norm(ki+1,1) = ki/kpts;
    %         ds.kx_norm(ki+1,1) = (((sqrt(3)/2)*ki/kpts)*(ki<kpts)+...                  % Gamma-X
    %                         (sqrt(3)/2)*(ki>=kpts && ki<2*kpts)+...              % X-M
    %                         (sqrt(3)/2)*(3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    %         ds.ky_norm(ki+1,1) = ((-1/2)*(ki<kpts)+...                          % Gamma-X
    %                         ((ki-kpts)/kpts-1/2)*(ki>=kpts && ki<2*kpts)+... % X-M
    %                         (1/2)*(3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    %     else
    %         ds.k_norm(ki+1,1) = ki/kpts;
    %         ds.kx_norm(ki+1,1) = ((ki/kpts)*(ki<kpts)+...                  % Gamma-X
    %                             1*(ki>=kpts && ki<2*kpts)+...              % X-M
    %                             (3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    %         ds.ky_norm(ki+1,1) = (0*(ki<kpts)+...                          % Gamma-X
    %                             (ki-kpts)/kpts*(ki>=kpts && ki<2*kpts)+... % X-M
    %                             (3*kpts-ki)/kpts*(ki>=2*kpts));            % M-Gamma
    %     end
    % else
    %     ds.k_norm = transpose(ds.kx/(pi/P.a));
    % end
    % save displacement and strain profile for all bands at 
    % high symmetry points (Gamma, X, M)
    if P.saveplots && mod(ki,kpts) == 0
        for nb = 1:nbands
            
            % text strings for k-points
            kpt_txts = {'G','X','M'};
            kpt_ttls = {'\Gamma','X','M'};
            
            % plot 3D displacement profile
            dispExpr = 'abs(solid.disp)';
            dispMax = mphmax(model,dispExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',ki+1);
            dispplot.set('looplevel',{num2str(nb) num2str(ki+1)});
            dispplot_vol.set('expr',dispExpr);
            dispplot_vol.set('rangedatamin','0').set('rangedatamax',dispMax);
            dispplot_vol.set('rangecolormin','0').set('rangecolormax',dispMax);
            dispplot.run;
            
            figure;
            mphplot(model,'dispplot');
            title(['Displacement, ',txt_sol,', ',kpt_ttls{ki/kpts+1},'-point, ',...
                   'band no. ',num2str(nb),', ',...
                   '\omega_m = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
            %title(['Displacement at ' kpt_txt '-point, band no. ' num2str(nb) ', eigenfreq = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
            colorbar('eastoutside');
            caxis([0 dispMax]);
            colormap(colortable('Rainbow'));
            daspect([1 1 1]);
            view(3)
            %axis tight
            axis off
            
            % filenames of displacement plots
            fname_fig = [txt_zsym,'_disp_',kpt_txts{ki/kpts+1},'_band_',num2str(nb),'.png'];
            pathFig = [P.datLoc,P.fileBase,'\',fname_fig];
            saveas(gcf,pathFig);
            close all;
        end  % for nb = 1:nbands
        
        % save fem data struct
        if ki == 0 
            ds.femG = fem;
        elseif ki == kpts 
            ds.femX = fem;
        elseif ki == 2*kpts %in: for nb = 1:nbands
            ds.femM = fem;
        end
        
    end  % P.saveplots && mod(ki,kpts) == 0
    
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
    ds.k_norm = ds.kx_norm;     
end

%% saving the mph files for debugging purposes 
if P.saveMPH
    path_mph = [P.datLoc,P.fileBase,'\bands_',txt_zsym,'.mph'];
    mphsave(model, path_mph);
end

% %% test save
% datLoc = 'L:\Individuals\cchia\OMC\COMSOL\FEMscripts\testCrossBands\';
% mphsave(model,[datLoc,'testCross.mph']);


ds.P = P;
end

