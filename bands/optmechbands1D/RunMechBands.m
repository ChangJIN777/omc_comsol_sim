% TFunction to compute mechanical bandstructure of a nanobeam unit cell.
% Modified from previous version (???) for compatibility with COMSOL 5.x.
%
% Input arguments: P
% P: data structure consisting of following parameters
%	Emod: Young's modulus; P.nu: Poisson's ratio; P.rho: density;
%   nbands: no. of bands to solve for; 
%   kpts: no. of k-points to iterate over (excluding gamma point);
%   freq: frequency at which to start looking for bands;
%   maxdof = max no. of DoFs in mesh;
%   eveny = 1 to solve for modes that are even in y;
%   saveplots = 1 to save displacement and strain profiles at G/X points
%   and other parameters for drawUnitCell_....m
%
% Output: ds 
% ds: data structure containing following variables
%   P: input parameters, with length of beam and rotated elasticity matrix
%   freqs: eigenfrequencies of nanobeam; 
%   F: (kpts+1) x nbands array of eigenfrequencies
%   kx: x-components of wavevectors
%   kx_norm: normalized x-components of wavevectors
%   femG, femX: Gamma- and X-point data structures containing 
%       equ: equation model with material data (rho, D, E, nu); 
%       bnds: boundary conditions (bcs), indices of boundaries
%             where periodic, symmetric, and anti-symmetric conditions are
%             applied (pbc_inds, sym_inds, asym_inds)
%       xmesh: extended mesh info
%       sol: eigenvalues (lambda, 1 x nbands array) and eigenfreqs 
%            (freq = -lambda/(2*pi*i), 1 x nbands array)
%
% Cleaven Chia, 09/29/16


function ds = RunMechBands(P,mfem,datLoc)

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);
ModelUtil.clear();
clear model
clear k_paramstr

% create COMSOL model named 'model' from which COMSOL methods can be called, 
% e.g. model.save
model = ModelUtil.create('model');

% utilFilesPath = 'G:\My Drive\OMC\COMSOL\FEMscripts20181224';
% addpath(utilFilesPath)

% constants
c = 299792458; % m/s
hbar = 1.05457148e-34; % J*s


% extract parameters from P
a = P.a;
wid = P.w;
thi = P.th;

if strcmp(P.xsect,'rect')
    thi = P.th/2;
%     evenz = P.mbevenz;
end

% if using triangular cross section, recalculate height of beam
if strcmp(P.xsect,'tri')
    thi = wid/(2*tan(P.theta*pi/180));
    P.mbevenz = 0;
end

nperiod = P.nperiod;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric



nbands = mfem.nbands;
kpts = mfem.kpts;
freq = mfem.mbfreq;
eveny = mfem.mbeveny;
evenz = mfem.mbevenz;

max_dof = mfem.max_dof;


% prefixes for filenames
if eveny == 1
    txt_ysym = 'eveny';
    txt_ysol = 'even y';
elseif eveny == -1
    txt_ysym = 'oddy';
    txt_ysol = 'odd y';
else
    txt_ysym = 'noysym';
    txt_ysol = 'no y sym';
end

if evenz == 1
    txt_zsym = 'evenz';
    txt_zsol = 'even z';
elseif evenz == -1
    txt_zsym = 'oddz';
    txt_zsol = 'odd z';
else
    txt_zsym = 'nozsym';
    txt_zsol = 'no z sym';
end


fBase = P.fileBase;

if ~exist([datLoc,fBase],'dir') && mfem.saveplots
    mkdir([datLoc,fBase])
end

geomname = 'ucell';

%% Define parameters for sweep over wavevectors
% parameter node
model.param.set('kx', '0');

% define k-points
k_liststr = ['0'];
for ki=0:kpts
    ds.kx(ki+1) = 0 + (pi/P.a)*ki/kpts;
    if ki > 0
        k_liststr = [k_liststr,', ',num2str(ds.kx(ki+1))];
    end
    k_paramstr{ki+1} = ['"k", "',num2str(ds.kx(ki+1)),'"'];
end

%% Set up geometry for unit cell
if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
    [model,P] = DrawBlockTet(model,P);
else
    [model,P] = DrawMechUnitCell(model,P,mfem);
end

if mfem.plotgeom
    figure;
    mphgeom(model)
    pathFig = [datLoc,fBase,'\',fBase,'_geom'];
    saveas(gcf,[pathFig,'.fig']);
    saveas(gcf,[pathFig,'.png']);
end
% beam = model.geom(geomname);    %handle for geometry node
% len = P.beamLen;
%display('Geometry done');

%% Define material and properties
P = LoadMaterialParams(P);

% E = P.E;
% nu = P.nu;
% rho = P.rho;

matTags = mphmodel(model.material);

% create material nodes
if ~isfield(matTags,P.beamMat)
    bMat = model.material.create(P.beamMat);
else
    bMat = model.material(P.beamMat);
end
bMat.label(P.beamMat);
bMat_def = bMat.propertyGroup('def');

% mechanical properties
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
mbfem.dia_domind = 1;%solid_index(beam,-1,[0,wid/2,thi],[0,0,1]);
bMat.selection.geom(geomname, 3);    %to select domain
bMat.selection.set(mbfem.dia_domind);

%% Set up physics (Solid Mechanics) and define boundary conditions
smech = model.physics.create('smech', 'SolidMechanics', 'ucell');
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
smech.selection.set(mbfem.dia_domind);
% lem = linear elastic material

% boundary conditions
% all BCs set to free by default
clear bnds

% periodic BC - set for faces in yz-plane at x = +/- a/2
perBCs = smech.create('perBCs', 'PeriodicCondition', 2);
perBCs.label('Periodic BC');
perBCs.set('PeriodicType', 'Floquet');
% bnds.pbc_inds = [];
% pbcinds_l = bndindex(beam, [0 0 0], [1 0 0]);
% pbcinds_r = bndindex(beam, [P.beamLen 0 0], [1 0 0]);
% pbc_inds = [pbcinds_l pbcinds_r];
% bnds.pbc_inds(end+1:end+length(pbc_inds)) = pbc_inds;
bnds.pbc_inds = [P.bndSel.Xsym,P.bndSel.Xend];
perBCs.selection.set(bnds.pbc_inds);
perBC_dest = perBCs.create('perBC_dest', 'DestinationDomains', 2);
perBC_dest.selection.set(bnds.pbc_inds(2));    %set face at x=a/2 as destination
perBCs.set('kFloquet', {'kx'; '0'; '0'});        %initialize Floquet vector

% symmetric BC
symBCs = smech.create('symBCs', 'SymmetrySolid', 2);
symBCs.label('Symmetric BC');
bnds.sym_inds = [];
if eveny == 1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.sym_inds = [bnds.sym_inds,P.bndSel.Ysym];
end
if evenz == 1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.sym_inds = [bnds.sym_inds,P.bndSel.Zsym];
end
if (~isempty(bnds.sym_inds))
    symBCs.selection.set(bnds.sym_inds);
end

% anti-symmetric BC
asymBCs = smech.create('asymBCs', 'Antisymmetry', 2);
asymBCs.label('Anti-symmetric BC');
bnds.asym_inds = [];
if eveny == -1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.asym_inds = [bnds.asym_inds,P.bndSel.Ysym];
end
if evenz == -1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.asym_inds = [bnds.asym_inds,P.bndSel.Zsym];
end
if (~isempty(bnds.asym_inds))
    asymBCs.selection.set(bnds.asym_inds);
end

mbfem.bnds = bnds;
disp('Solid Mechanics added - boundary conditions done');

%% Add study and solver sequences
% add parametric and eigenfrequency study
study = model.study.create('study');
std_param = study.create('std_param', 'Parametric');
std_param.set('pname', 'kx');
std_param.set('plistarr', k_liststr);
std_param.set('punit', '');
std_eigv = study.create('std_eigv','Eigenfrequency');
std_eigv.set('neigsactive',true).set('neigs',nbands);
std_eigv.set('shiftactive',true).set('shift',num2str(0-1i*2*pi*freq));

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
solv_eigv.set('eigref',num2str(0-1i*2*pi*freq));
solv_eigv.feature('dDef').set('linsolver', 'spooles');
solv_eigv.feature('aDef').set('complexfun', 'off');

psolv = model.sol.create('psolv');
psolv.study('study');

% add batch job configuration for parameter sweep
pbatch = model.batch.create('pbatch', 'Parametric');
pbatch_solseq = pbatch.create('pbatch_solseq', 'Solutionseq');
pbatch.study('study');
pbatch.attach('study');
pbatch.set('pname', 'kx');
pbatch.set('plistarr', k_liststr);
pbatch.set('punit', '');
pbatch.set('err', true);
pbatch.set('control', 'std_param');
pbatch_solseq.set('psol', 'psolv');
pbatch_solseq.set('param', k_paramstr);
pbatch_solseq.set('seq', 'solv');

disp('Study, solver, sweep and batch nodes added');

%% Add Mesh
fem = mbfem;
fem.mesh_quality = mfem.meshSize;
mesh = model.mesh.create('mesh', 'ucell');
display(['Meshing with quality: ' num2str(fem.mesh_quality)]);
mesh.autoMeshSize(fem.mesh_quality).run;
%mphmesh(beammodel);    %to plot mesh

%% Adjust mesh if max DOFs exceeded
mesh_ok = 0;
% Xmesh info - get degrees of freedom
fem.xmesh = mphxmeshinfo(model, 'soltag', 'solv', ...
                                   'studysteptag', 'solv_stdstep');
dofs = fem.xmesh.ndofs;
display(['Estimated no. of DoFs: ' num2str(dofs)]);

while (~mesh_ok) && (fem.mesh_quality < 10)
    while (dofs >= max_dof) && (fem.mesh_quality < 10)
        fem.mesh_quality = fem.mesh_quality + 1;
        display(['Meshing with quality: ' num2str(fem.mesh_quality)]);
        mesh.autoMeshSize(fem.mesh_quality).run;
        fem.xmesh = mphxmeshinfo(model, 'soltag', 'solv', ...
                                           'studysteptag', 'solv_stdstep');
        dofs = fem.xmesh.ndofs;
        display(['Estimated no. of DoFs: ' num2str(dofs)]);
    end
    mesh_ok = 1;
end

%% Solve for bands
solv.runAll;
pbatch.run;
%display('Mechanical eigenmodes solved');

%% Set up results node
% change name of default dataset
model.result.dataset('dset1').tag('dset');
dset = model.result.dataset('dset');
dset.set('solution', 'solv');
model.result.dataset('dset2').tag('pdset');
pdset = model.result.dataset('pdset');
pdset.set('solution', 'psolv');

if mfem.saveplots
    % create sector dataset
    xySymFac = 2^abs(eveny);
    zSym = abs(evenz)*strcmp(P.xsect,'rect');
    zSymFac = 2^zSym;
    
    dsetTags = mphmodel(model.result.dataset);
    
    if ~isfield(dsetTags,'pdset_sec')
        pdset_sec = model.result.dataset.create('pdset_sec', 'Sector3D');
    else
        pdset_sec = model.result.dataset('pdset_sec');
    end
    pdset_sec.label('M Sol Full Beam XY');
    pdset_sec.set('data', 'pdset');
    pdset_sec.set('method', 'twopoint');
    pdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
    pdset_sec.setIndex('genpoints', '0', 0, 1);
    pdset_sec.setIndex('genpoints', '0', 0, 2);
    pdset_sec.setIndex('genpoints', '0', 1, 0);
    pdset_sec.setIndex('genpoints', '0', 1, 1);
    pdset_sec.setIndex('genpoints', '1', 1, 2);
    pdset_sec.set('sectors', xySymFac);
    pdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
    pdset_sec.set('reflaxis', {'1' '0' '0'});   % reflection axis
    if eveny==-1
        pdset_sec.set('rotinv', 'on');  % odd symmetry about x-plane
        pdset_sec.set('reflinv', 'on');
    end
    pdsetPlot = 'pdset_sec';
    
    if zSym
        if ~isfield(dsetTags,'pdset_secZ')
            pdset_sec = model.result.dataset.create('pdset_secZ', 'Sector3D');
        else
            pdset_sec = model.result.dataset('pdset_secZ');
        end
        pdset_sec.label('M Sol Full Beam XYZ');
        pdset_sec.set('data', 'pdset_sec');
        pdset_sec.set('method', 'twopoint');
        pdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
        pdset_sec.setIndex('genpoints', '0', 0, 1);
        pdset_sec.setIndex('genpoints', '0', 0, 2);
        pdset_sec.setIndex('genpoints', '1', 1, 0);
        pdset_sec.setIndex('genpoints', '0', 1, 1);
        pdset_sec.setIndex('genpoints', '0', 1, 2);
        pdset_sec.set('sectors', zSymFac);
        pdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
        pdset_sec.set('reflaxis', {'0' '1' '0'});   % reflection axis
        if evenz==-1
            pdset_sec.set('rotinv', 'on');  % odd symmetry about z-plane
            pdset_sec.set('reflinv', 'on');
        end
        pdsetPlot = 'pdset_secZ';
        
    end
    
%     model.result.dataset.create('sec1', 'Sector3D');
%     model.result.dataset('sec1').set('trans', 'rotrefl');
%     model.result.dataset('sec1').set('reflaxis', {'1' '0' '0'});
%     model.result.dataset('sec1').set('data', 'pdset');
%     if ~eveny
%         model.result.dataset('sec1').set('reflinv', 'on');
%     end
    
    
    % create 3D plot group for displacement field
    dispplot = model.result.create('dispplot', 'PlotGroup3D');
    dispplot.set('data',pdsetPlot);
    dispplot.set('solrepresentation', 'solutioninfo');
    dispplot.set('titletype', 'none');
    dispplot_vol = dispplot.create('dispplot_vol', 'Volume');
    dispplot_vol.set('rangedataactive','on').set('rangecoloractive','on');
    dispplot_def = dispplot_vol.create('dispplot_def', 'Deform');
    dispplot_vol.set('data', 'parent');
    dispplot.run;

    % create 3D plot group for strain profile
    strplot = model.result.create('strplot', 'PlotGroup3D');
    strplot.set('data',pdsetPlot);
    strplot.set('solrepresentation', 'solutioninfo');
    strplot.set('titletype', 'none');
    strplot_slc = strplot.create('strplot_slc', 'Slice');
    strplot_slc.set('planetype', 'quick');
    strplot_slc.set('quickxmethod', 'coord').set('quickx', (1-holeatedge/2)*a);
    strplot_slc.set('rangedataactive','on').set('rangecoloractive','on');
    strplot_slc.set('data', 'parent');
    strplot.run;
    
    if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
        strplot2 = model.result.create('strplot2', 'PlotGroup3D');
        strplot2.set('data',pdsetPlot);
        strplot2.set('solrepresentation', 'solutioninfo');
        strplot2.set('titletype', 'none');
        strplot2_slc = strplot2.create('strplot2_slc', 'Slice');
        strplot2_slc.set('planetype', 'quick');
        strplot2_slc.set('quickxmethod', 'coord').set('quickx', (1-(1-holeatedge)/2)*a);
        strplot2_slc.set('rangedataactive','on').set('rangecoloractive','on');
        strplot2_slc.set('data', 'parent');
        strplot2.run;
    end
    
end

%% Save data and plots

% extract solution info from parameter sweep
sols = mphsolutioninfo(model);
lambda_inds = find(strcmp(sols.psolv.mapheaders,'lambda'));
kx_inds = find(strcmp(sols.psolv.mapheaders,'kx'));
inner_inds = find(strcmp(sols.psolv.mapheaders,'Inner'));
outer_inds = find(strcmp(sols.psolv.mapheaders,'Outer'));

% assemble solutions
for ki = 0:kpts
    % assemble eigenvalues and eigenfrequencies
    lambda_ki = find(sols.psolv.map(:,outer_inds)==ki+1);
    fem.sol.lambda = sols.psolv.map(lambda_ki,lambda_inds);
    fem.sol.freqs = abs(fem.sol.lambda)/(2*pi);
    for nb = 1:nbands
        ds.F(ki+1,nb) = fem.sol.freqs(nb);
    end
    ds.kx_norm = transpose(ds.kx/(pi/P.a));
    
    % assenmble eigenvector array - needed? kiv
    
    % save displacement and strain profile for all bands at gamma or X points
    if mfem.saveplots && (ki==0 || ki==kpts)   %in: for ki = 0:kpts
        for nb = 1:nbands %in: if P.saveplots && (ki==0 || ki==kpts)
            
            % strings for k-points
            if ki == 0  %in: for nb = 1:nbands
                kpt_txt = 'G';
            else
                kpt_txt = 'X';
            end
            
            dispExpr = 'abs(solid.disp)';
            dispMax = mphmax(model,dispExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',ki+1);
            ds.dispMax(ki+1,nb) = dispMax;
                         
            % save strains for gamma-point, symmetric array
%             if ki == 0 && eveny == 1 && evenz == 1
%                 delta = 10e-9;
%                 xSiV = linspace(0,P.a,ceil(P.a/delta));
%                 ySiV = linspace(0,P.w/2,ceil(P.w/delta));
%                 zSiV = linspace(0,P.th/2,ceil(P.th/delta));
%                 
%                 [xm,ym,zm] = meshgrid(xSiV,ySiV,zSiV); % array of size length(yCoords) x length(xCoords(xi)) x length(zCoords)
%                 coord = transpose([xm(:),ym(:),zm(:)]);   % array of size 3x(length(xCoords(xi))*length(yCoords)*length(zCoords))
%                 
%                 [eXX,eYY,eZZ,eXY,eYZ,eXZ] = mphinterp(model,...
%                 {'solid.eXX','solid.eYY','solid.eZZ', ...
%                     'solid.eXY','solid.eYZ','solid.eXZ'},'coord',coord,...
%                     'dataset','pdset','edim',3,'selection','all','solnum',nb,'outersolnum',ki+1);
%                 for ci = 1:size(coord,2)
%                     eTens = [eXX(ci), eXY(ci), eXZ(ci);
%                             eXY(ci), eYY(ci), eYZ(ci);
%                             eXZ(ci), eYZ(ci), eZZ(ci)];
% 
%                     % then rotate into [111] and [-111] bases
%                     [e111, em111] = rotateStrainSiV(eTens,P.rxtal);
% 
%                     % extract relevant strain components from strain tensors in [111] and [-111] bases
%                     e111YY = e111(2,2); 
%                     e111ZZ = e111(3,3);
%                     e111XY = e111(1,2); 
%                     e111XZ = e111(1,3);
%                     e111YZ = e111(2,3);
% 
%                     em111YY = em111(2,2);
%                     em111ZZ = em111(3,3);
%                     em111XY = em111(1,2);
%                     em111XZ = em111(1,3);
%                     em111YZ = em111(2,3);
% 
%                     % calc transverse strain coupling (not normalized)
%                     betaG = P.dg*(e111YY-e111ZZ) + P.fg*e111XY;
%                     gammaG = -2*P.dg*e111YZ + P.fg*e111XZ;
%                     lambdaG_111(nb,ci) = sqrt(real(betaG)^2+real(gammaG)^2);
% 
%                     betaGm = P.dg*(em111YY-em111ZZ) + P.fg*em111XY;
%                     gammaGm = -2*P.dg*em111YZ + P.fg*em111XZ;
%                     lambdaG_m111(nb,ci) = sqrt(real(betaGm)^2+real(gammaGm)^2);
%                 end
%                 
%                 % find max for each mode
%                 [lambdaG_111Max, lambdaG_111MaxIdx] = max(lambdaG_111,[],2);
%                 [lambdaG_m111Max, lambdaG_m111MaxIdx] = max(lambdaG_m111,[],2);
%                 
%                 lambdaG_111MaxCoord = coord(:,lambdaG_111MaxIdx);
%                 lambdaG_m111MaxCoord = coord(:,lambdaG_m111MaxIdx);
%                 
%                 % save to data struct
%                 ds.lambdaG_111Max = lambdaG_111Max;
%                 ds.lambdaG_111MaxIdx = lambdaG_111MaxIdx;
%                 ds.lambdaG_111MaxCoord = lambdaG_111MaxCoord;
%                 
%                 ds.lambdaG_m111Max = lambdaG_m111Max;
%                 ds.lambdaG_m111MaxIdx = lambdaG_m111MaxIdx;
%                 ds.lambdaG_m111MaxCoord = lambdaG_m111MaxCoord;
%                 
%                 % plot
% %                 figure; mphgeom(model,'ucell','facealpha',0); box on
%                 
%             end
            
            % plot 3D displacement profile
            
            dispplot.set('looplevel',{num2str(nb) num2str(ki+1)});
            dispplot_vol.set('expr',dispExpr);
            dispplot_vol.set('rangedatamin','0').set('rangedatamax',dispMax);
            dispplot_vol.set('rangecolormin','0').set('rangecolormax',dispMax);
            dispplot.run;
            
            figure;
            try %in: for nb = 1:nbands
                mphplot(model,'dispplot');
            catch err
            end
            title(['Displacement, ',txt_ysol,', ',txt_zsol,', ',kpt_txt,'-point, ',...
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
            fname_fig = [txt_ysym,'_',txt_zsym,'_disp_',kpt_txt,'_band_',num2str(nb),'.png'];
            pathFig = [datLoc,fBase,'\',fname_fig];
            saveas(gcf,pathFig);
                        
            % plot cross-sectional strain profile
            strExpr = 'ucell.solid.eXX+ucell.solid.eYY+ucell.solid.eZZ';
            absStrExpr = ['abs(' strExpr ')'];
            strMax = mphmax(model,absStrExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',ki+1);
            ds.strMax(ki+1,nb) = strMax;
            %strplot.set('solnum', nb);
            strplot.set('looplevel',{num2str(nb) num2str(ki+1)});
            strplot_slc.set('expr', strExpr);
            strplot_slc.set('rangedatamin', -strMax).set('rangedatamax', strMax);
            strplot_slc.set('rangecolormin', -strMax).set('rangecolormax', strMax);
            strplot_slc.set('colortable', 'WaveLight');
            strplot.run;
            
            figure;
            try
                mphplot(model,'strplot');
            catch err
            end
            title(['Strain, ',txt_ysol,', ',txt_zsol,', ',kpt_txt,'-point, ',...
                   'band no. ',num2str(nb),', ',...
                   '\omega_m = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
            %title(['Strain at ' kpt_txt '-point, band no. ' num2str(nb) ', eigenfreq = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
            colorbar('eastoutside');
            caxis([-strMax strMax]);
            colormap(colortable('WaveLight'));
            daspect([1 1 1]);
            view(90,0)
            %axis tight
            axis off
            
            % filenames of strain plots
            fname_fig = [txt_ysym,'_',txt_zsym,'_strain_',kpt_txt,'_band_',num2str(nb),'.png'];
            pathFig = [datLoc,fBase,'\',fname_fig];
            saveas(gcf,pathFig);
            
            if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
                strplot2.set('looplevel',{num2str(nb) num2str(ki+1)});
                strplot2_slc.set('expr', strExpr);
                strplot2_slc.set('rangedatamin', -strMax).set('rangedatamax', strMax);
                strplot2_slc.set('rangecolormin', -strMax).set('rangecolormax', strMax);
                strplot2_slc.set('colortable', 'WaveLight');
                strplot.run;

                figure;
                try
                    mphplot(model,'strplot2');
                catch err
                end
                title(['Strain, ',txt_ysol,', ',txt_zsol,', ',kpt_txt,'-point, ',...
                       'band no. ',num2str(nb),', ',...
                       '\omega_m = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
                %title(['Strain at ' kpt_txt '-point, band no. ' num2str(nb) ', eigenfreq = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
                colorbar('eastoutside');
                caxis([-strMax strMax]);
                colormap(colortable('WaveLight'));
                daspect([1 1 1]);
                view(90,0)
                %axis tight
                axis off

                % filenames of strain plots
                fname_fig = [txt_ysym,'_',txt_zsym,'_strain2_',kpt_txt,'_band_',num2str(nb),'.png'];
                pathFig = [datLoc,fBase,'\',fname_fig];
                saveas(gcf,pathFig);
            end
            
            close all;
        end %of: for nb = 1:nbands
        
        % save fem data struct
        if ki == 0 %in: for nb = 1:nbands
            ds.femG = fem;
        end
        
        if ki == kpts %in: for nb = 1:nbands
            ds.femX = fem;
        end


    end %of: if P.saveplots && (ki==0 || ki==kpts) 
    
    

end %of: for ki = 0:kpts
%display('Postprocessing done');

if mfem.savedat
    mphfname = ['Bands_',txt_ysym,'_',txt_zsym,'.mph'];
    path_mph = [datLoc,fBase,'\',mphfname];
    mphsave(model, path_mph);
end
ModelUtil.clear();
%display('Model and data saved. Run complete');
% rmpath(utilFilesPath)
end
