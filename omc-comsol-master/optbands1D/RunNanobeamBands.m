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


function ds = RunNanobeamBands(P,datLoc)

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);
ModelUtil.clear();
clear model

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
    P.obevenz = 0;
end

nperiod = P.nperiod;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric
P.nholes = nperiod + holeatedge;


nbands = P.nbands;
kpts = P.kpts;
freq = P.obfreq;
eveny = P.obeveny;
evenz = P.obevenz;

max_dof = P.max_dof;


% prefixes for filenames
if strcmp(P.xsect,'tri') || (strcmp(P.xsect,'rect') && ~abs(evenz))
    if eveny == 1
        txt_sym = 'TM';
        txt_sol = 'TM mode';
    elseif eveny == -1
        txt_sym = 'TE';
        txt_sol = 'TE mode';
    end
elseif strcmp(P.xsect,'rect') && abs(evenz)
    if eveny == 1 && evenz == -1
        txt_sym = 'TM';
        txt_sol = 'TM mode';
    elseif eveny == -1 && evenz == 1
        txt_sym = 'TE';
        txt_sol = 'TE mode';
    end
end
% if eveny == 1
%     txt_ysym = 'eveny';
%     txt_ysol = 'even y';
% elseif eveny == -1
%     txt_ysym = 'oddy';
%     txt_ysol = 'odd y';
% else
%     txt_ysym = 'noysym';
%     txt_ysol = 'no y sym';
% end
% 
% if evenz == 1
%     txt_zsym = 'evenz';
%     txt_zsol = 'even z';
% elseif evenz == -1
%     txt_zsym = 'oddz';
%     txt_zsol = 'odd z';
% else
%     txt_zsym = 'nozsym';
%     txt_zsol = 'no z sym';
% end

% create base folder name
if ~isfield(P,'fileBase')
    fBase = ['a',num2str(P.a*1e9,'%.0f'),'nm_',...
        'w',num2str(P.w*1e9,'%.0f'),'nm_', ...
        'hx',num2str(P.hx*1e9,'%.0f'),'nm_', ...
        'hy',num2str(P.hy*1e9,'%.0f'),'nm_',...
        'rxt',num2str(P.rxtal),'o'];
    if strcmp(P.xsect,'tri')
        fBase = [num2str(P.theta,'%.0f'),'o_',fBase];
    elseif strcmp(P.xsect,'rect')
        fBase = ['th',num2str(P.th*1e9,'%.0f'),'nm_',fBase];
    end
    if isfield(P,'prefname')
        fBase = [P.prefname,'_',fBase];
    end
    
    P.fileBase = fBase;
end
fBase = P.fileBase;

if ~exist([datLoc,fBase],'dir') && P.saveplots
    mkdir([datLoc,fBase])
end

geomname = 'ucell';

%% Define parameters for sweep over wavevectors
% parameter node
model.param.set('kx', '0');

% define k-points
k_liststr = [];
k_start = 20;
k_int = 1;


k_idxs = k_start:k_int:kpts;
ds.kx = (pi/P.a)*k_idxs/kpts;
for kn = 1:length(k_idxs)
    ki = k_idxs(kn);
    if isempty(k_liststr)
        k_liststr = num2str(ds.kx(kn));
    else
        k_liststr = [k_liststr,', ',num2str(ds.kx(kn))];
    end
    k_paramstr{kn} = ['"k", "',num2str(ds.kx(kn)),'"'];
end

%% Set up geometry for unit cell
if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
    [model,P] = DrawBlockTet(model,P);
else
    [model,P] = DrawUnitCell(model,P);
end

if P.plotgeom
    figure;
    mphgeom(model,'ucell','facealpha',0.5)
    pathFig = [datLoc,fBase,'\',fBase,'_geom'];
    saveas(gcf,[pathFig,'.fig']);
    saveas(gcf,[pathFig,'.png']);
end
% beam = model.geom(geomname);    %handle for geometry node
% len = P.beamLen;
%display('Geometry done');

%% Define material and properties
P = LoadMaterialParams(P);
nbeam = P.nbeam;

matTags = mphmodel(model.material);

% create material nodes
if ~isfield(matTags,P.beamMat)
    bMat = model.material.create(P.beamMat);
else
    bMat = model.material(P.beamMat);
end
bMat.label(P.beamMat);
bMat_def = bMat.propertyGroup('def');

% optical properties - conductivity, rel permittivity, refractive index
% first/second element of cell array refers to value for air/beam (e.g. diamond)
% create air material node, if not already created
if ~isfield(matTags,'air')
    air = model.material.create('air');
else
    air = model.material('air');
end
air.label('Air');
air_def = air.propertyGroup('def');

obfem.sigma = {0,'1e-12[S/m]'};  % conductivity
obfem.epsilonr = {1,nbeam^2};    % relative permittivity
obfem.n = {1,nbeam};             % refractive index

% define optical properties of air and beam materials
air_ref = air.propertyGroup.create('air_ref', 'Refractive index');
air_def.set('electricconductivity', obfem.sigma{1});
air_def.set('relpermittivity', obfem.epsilonr{1});
air_def.set('relpermeability', 1);
air_ref.set('n', obfem.n{1});

bMat_ref = bMat.propertyGroup.create('bMat_ref', 'Refractive index');
bMat_def.set('electricconductivity', obfem.sigma{2});
bMat_def.set('relpermittivity', obfem.epsilonr{2});
bMat_def.set('relpermeability', 1);
bMat_ref.set('n', obfem.n{2});

% apply material to each domain
obfem.dia_domind = P.domSel.ucell;
bMat.selection.geom(geomname, 3);    %to select domain
bMat.selection.set(obfem.dia_domind);

obfem.air_domind = P.domSel.cyl;
air.selection.geom(geomname, 3);    %to select domain
air.selection.set(obfem.air_domind);

display(['Beam material: ',P.beamMat,' added']);

%% Set up physics (RF EMW freq domain) and define boundary conditions
emw = model.physics.create('emw', 'ElectromagneticWaves', geomname);
emw.selection.set([obfem.dia_domind, obfem.air_domind]);

% boundary conditions: default = PEC for boundaries adjacent to single domain
% - not applicable for boundaries adjacent to two domains (i.e. between nanobeam and air cylinder)
% (this is the continuity condition in COMSOL 3.5)
clear bnds

% periodic BC - set for faces in yz-plane at x = +/- a/2
perBCs = emw.create('perBCs', 'PeriodicCondition', 2);
perBCs.label('Periodic BC');
perBCs.set('PeriodicType', 'Floquet');
bnds.pbc_inds = [P.bndSel.Xsym,P.bndSel.Xend];
perBCs.selection.set(bnds.pbc_inds);
perBC_dest = perBCs.create('perBC_dest', 'DestinationDomains', 2);
perBC_dest.selection.set(P.bndSel.Xend);    %set face at x=a/2 as destination
perBCs.set('kFloquet', {'kx'; '0'; '0'});        %initialize Floquet vector

% Perfect magnetic conductor (PMC) - set for even modes
PMCbnd = emw.create('PMCbnd', 'PerfectMagneticConductor', 2);
PMCbnd.label('Perfect Magnetic Conductor');
bnds.pmc_inds = [];
if eveny == 1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.pmc_inds = [bnds.pmc_inds,P.bndSel.Ysym];
end
if evenz == 1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.pmc_inds = [bnds.pmc_inds,P.bndSel.Zsym];
end
if ~isempty(bnds.pmc_inds)
    PMCbnd.selection.set(bnds.pmc_inds);
else
    PMCbnd.active(false);
end


% Scattering boundary condition (SBC) - apply to ends of beam and curved bnd
% of air cylinder
scatbnd = emw.create('scatbnd', 'Scattering', 2);
scatbnd.label('Scattering Boundary Condition');
bnds.scat_inds = P.bndSel.cylCurv;
if ~isempty(bnds.scat_inds)
    scatbnd.selection.set(bnds.scat_inds);
else
    scatbnd.active(false);
end

% Perfect electric conductor (PEC) - set for odd modes
% compute boundary indices to check if these are the same as the 
% remaining PEC boundaries set by default
emw.feature('pec1').tag('PECbnd');
PECbnd = emw.feature('PECbnd');
PECbnd.label('Perfect Electric Conductor');
bnds.pec_inds = [];
if eveny == -1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.pec_inds = [bnds.pec_inds,P.bndSel.Ysym];
end
if evenz == -1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.pec_inds = [bnds.pec_inds,P.bndSel.Zsym];
end
if ~isempty(bnds.pec_inds)
    pec_inds_default = transpose(PECbnd.selection.entities(2));
    if ~isempty(setdiff(bnds.pec_inds, pec_inds_default))
        %display(['bnds from bndindex:' num2str(bnds.pec_inds)]);
        %display(['bnds from default:' num2str(pec_inds_default)]);
        error('perfect electric conductor - boundaries not set correctly');
%     else
%         display('boundaries set correctly');
%         display(['perfect electric conductor - bnds:' num2str(bnds.pec_inds)]);
    end
end

obfem.bnds = bnds;

display('EM Waves Frequency Domain added - boundary conditions done');

%% Add study and solver sequences
% add parametric and eigenfrequency study
study = model.study.create('study');
std_param = study.create('std_param', 'Parametric');
std_param.set('pname', 'kx');
std_param.set('plistarr', k_liststr);
std_param.set('punit', '');
std_eigv = study.create('std_eigv','Eigenfrequency');
std_eigv.set('neigsactive',true).set('neigs',nbands).set('eigunit', 'Hz');
std_eigv.set('shiftactive',true).set('shift',num2str(freq));
std_eigv.set('eigwhich', 'lr'); % find eigenvalues with larger real part compared to shift

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
pbatch.set('pname', 'kx');
pbatch.set('plistarr', k_liststr);
pbatch.set('punit', '');
pbatch.set('err', true);
pbatch.set('control', 'std_param');
pbatch_solseq.set('psol', 'psolv');
pbatch_solseq.set('param', k_paramstr);
pbatch_solseq.set('seq', 'solv');

display('Study, solver, sweep and batch nodes added');

%% Add Mesh
fem = obfem;
fem.mesh_quality = 5;
mesh = model.mesh.create('mesh', 'ucell');
mesh.feature('size').set('hauto', fem.mesh_quality);
display(['Meshing with quality: ' num2str(fem.mesh_quality)]);
% mesh.autoMeshSize(fem.mesh_quality).run;

mesh.create('ftri2', 'FreeTri');
mesh.feature('ftri2').selection.set([P.bndSel.Xsym]);
mesh.run;

mesh.create('cpf1', 'CopyFace');
mesh.feature('cpf1').selection('source').set(P.bndSel.Xsym);
mesh.feature('cpf1').selection('destination').set(P.bndSel.Xend);
mesh.run;

mesh.create('ftet1', 'FreeTet');
mesh.feature('ftet1').selection.geom('ucell', 3);
mesh.feature('ftet1').selection.set([1 2]);
mesh.feature('ftet1').create('size1', 'Size');
mesh.feature('ftet1').feature('size1').set('hauto', fem.mesh_quality);
mesh.run;
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

if P.savedat
    mphfname = ['Bands_',txt_sym,'.mph'];
    path_mph = [datLoc,fBase,'\',mphfname];
    mphsave(model, path_mph);
end

%% Solve for bands
solv.runAll;
pbatch.run;

%% Set up results node
% change name of default dataset
model.result.dataset('dset1').tag('dset');
dset = model.result.dataset('dset');
dset.set('solution', 'solv');
model.result.dataset('dset2').tag('pdset');
pdset = model.result.dataset('pdset');
pdset.set('solution', 'psolv');

%%
if P.saveplots
    % create sector dataset for each parity per symmetry plane
    xySymFac = 2^abs(eveny);
    zSym = abs(evenz)*strcmp(P.xsect,'rect');
    zSymFac = 2^zSym;
    
    dsetTags = mphmodel(model.result.dataset);
    
    % even Y
    if ~isfield(dsetTags,'sec_Yeven')
        secdset = model.result.dataset.create('sec_Yeven', 'Sector3D');
    else
        secdset = model.result.dataset('sec_Yeven');
    end
    secdset.label('Full unit cell - Y even');
    secdset.set('data', 'pdset');
    secdset.set('method', 'twopoint');
    secdset.setIndex('genpoints', '0', 0, 0); % set axis data
    secdset.setIndex('genpoints', '0', 0, 1);
    secdset.setIndex('genpoints', '0', 0, 2);
    secdset.setIndex('genpoints', '0', 1, 0);
    secdset.setIndex('genpoints', '0', 1, 1);
    secdset.setIndex('genpoints', '1', 1, 2);
    secdset.set('sectors', xySymFac);
    secdset.set('trans', 'rotrefl');          % transformation: rotate and reflect
    secdset.set('reflaxis', {'1' '0' '0'});   % reflection axis
    secdset.set('rotinv', 'on');
    secdset.set('reflinv', 'off');
    
    % odd Y
    if ~isfield(dsetTags,'sec_Yodd')
        secdset = model.result.dataset.duplicate('sec_Yodd','sec_Yeven');
    else
        secdset = model.result.dataset('sec_Yodd');
    end
    secdset.label('Full unit cell - Y odd');
    secdset.set('reflinv', 'on');
    
    % designate datasets for each plot
    
    if ~zSym
        EnormSec = 'sec_Yeven';
        if eveny == 1
            EySec = 'sec_Yodd';
            EzSec = 'sec_Yeven';
        elseif eveny == -1
            EySec = 'sec_Yeven';
            EzSec = 'sec_Yodd';
        end
    elseif zSym
        % even Y odd Z
        if ~isfield(dsetTags,'sec_YevenZodd')
            secdset = model.result.dataset.duplicate('sec_YevenZodd','sec_Yodd');
        else
            secdset = model.result.dataset('sec_YevenZodd');
        end
        secdset.label('Full unit cell - Y even Z odd');
        secdset.set('data', 'sec_Yeven');
        secdset.setIndex('genpoints', '1', 1, 0);
        secdset.setIndex('genpoints', '0', 1, 1);
        secdset.setIndex('genpoints', '0', 1, 2);
        secdset.set('sectors', zSymFac);
        secdset.set('trans', 'rotrefl');          % transformation: rotate and reflect
        secdset.set('reflaxis', {'0' '1' '0'});   % reflection axis
        secdset.set('reflinv', 'off');
        
        % odd Y even Z
        if ~isfield(dsetTags,'sec_YoddZeven')
            secdset = model.result.dataset.duplicate('sec_YoddZeven','sec_Yodd');
        else
            secdset = model.result.dataset('sec_YoddZeven');
        end
        secdset.label('Full unit cell - Y odd Z even');
        secdset.set('data', 'sec_Yodd');
        secdset.setIndex('genpoints', '1', 1, 0);
        secdset.setIndex('genpoints', '0', 1, 1);
        secdset.setIndex('genpoints', '0', 1, 2);
        secdset.set('sectors', zSymFac);
        secdset.set('trans', 'rotrefl');          % transformation: rotate and reflect
        secdset.set('reflaxis', {'0' '1' '0'});   % reflection axis
        secdset.set('reflinv', 'on');
        
        % designate datasets for each plot
        EnormSec = 'sec_YevenZodd';
        if eveny == 1 && evenz == -1
            EySec = 'sec_YoddZeven';
            EzSec = 'sec_YevenZodd';
        elseif eveny == -1 && evenz == 1
            EySec = 'sec_YevenZodd';
            EzSec = 'sec_YoddZeven';
        end
    end
    
    % create plot group for E-field intensity
    Enormplt = model.result.create('Enormplt', 'PlotGroup3D');
    Enormplt.set('data',EnormSec);
    Enormplt.set('solrepresentation', 'solutioninfo');
    Enormplt.set('titletype', 'none');
    Enormplt_slc = Enormplt.create('Enormplt_slc', 'Slice');
    Enormplt_slc.set('planetype', 'quick');
    Enormplt_slc.set('quickxmethod', 'coord').set('quickx', (1-holeatedge/2)*a);
    Enormplt_slc.set('rangedataactive','on').set('rangecoloractive','on');
    Enormplt_slc.set('data', 'parent');
    Enormplt.run;

    % create 3D plot group for Ey
    Eyplt = model.result.create('Eyplt', 'PlotGroup3D');
    Eyplt.set('data',EySec);
    Eyplt.set('solrepresentation', 'solutioninfo');
    Eyplt.set('titletype', 'none');
    Eyplt_slc = Eyplt.create('Eyplt_slc', 'Slice');
    Eyplt_slc.set('planetype', 'quick');
    Eyplt_slc.set('quickxmethod', 'coord').set('quickx', (1-holeatedge/2)*a);
    Eyplt_slc.set('rangedataactive','on').set('rangecoloractive','on');
    Eyplt_slc.set('data', 'parent');
    Eyplt.run;
    
    % create 3D plot group for Ey
    Ezplt = model.result.create('Ezplt', 'PlotGroup3D');
    Ezplt.set('data',EzSec);
    Ezplt.set('solrepresentation', 'solutioninfo');
    Ezplt.set('titletype', 'none');
    Ezplt_slc = Ezplt.create('Ezplt_slc', 'Slice');
    Ezplt_slc.set('planetype', 'quick');
    Ezplt_slc.set('quickxmethod', 'coord').set('quickx', (1-holeatedge/2)*a);
    Ezplt_slc.set('rangedataactive','on').set('rangecoloractive','on');
    Ezplt_slc.set('data', 'parent');
    Ezplt.run;
    
    
end

%% Save data and plots

% extract solution info from parameter sweep
sols = mphsolutioninfo(model);
lambda_inds = find(strcmp(sols.psolv.mapheaders,'lambda'));
kx_inds = find(strcmp(sols.psolv.mapheaders,'kx'));
inner_inds = find(strcmp(sols.psolv.mapheaders,'Inner'));
outer_inds = find(strcmp(sols.psolv.mapheaders,'Outer'));

%% assemble solutions
ds.kx_norm = transpose(ds.kx/(pi/P.a));
for kn = 1:length(k_idxs)
    ki = k_idxs(kn);
    
    % assemble eigenvalues and eigenfrequencies
    lambda_ki = find(sols.psolv.map(:,outer_inds)==kn);
    fem.sol.lambda = sols.psolv.map(lambda_ki,lambda_inds);
    fem.sol.freqs = abs(fem.sol.lambda)/(2*pi);
    for nb = 1:nbands
        ds.F(kn,nb) = fem.sol.freqs(nb);
    end
     
    % save displacement and strain profile for all bands at gamma or X points
    if P.saveplots && ki==kpts   
        for nb = 1:nbands 
            
            % strings for k-points
            kpt_txt = 'X';
            
            EnormExpr = 'emw.normE';
            EnormMax = mphmax(model,EnormExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',kn);
            ds.EnormMax(kn,nb) = EnormMax;

            % plot 3D displacement profile
            
            Enormplt.set('looplevel',{num2str(nb) num2str(kn)});
            Enormplt_slc.set('expr',EnormExpr);
            Enormplt_slc.set('rangedatamin','0').set('rangedatamax',EnormMax);
            Enormplt_slc.set('rangecolormin','0').set('rangecolormax',EnormMax);
            Enormplt.run;
            
            figure;
            try %in: for nb = 1:nbands
                mphplot(model,'Enormplt');
            catch err
            end
            title({['normalized E-field, ',txt_sol,', ',kpt_txt,'-point, '];...
                   ['band no. ',num2str(nb),', ',...
                   '\omega_o = ' num2str(ds.F(kn,nb)/1e12) ' THz']});
            %title(['Displacement at ' kpt_txt '-point, band no. ' num2str(nb) ', eigenfreq = ' num2str(ds.F(ki+1,nb)/1e9) ' GHz']);
            colorbar('eastoutside');
            caxis([0 EnormMax]);
            colormap(colortable('Rainbow'));
            daspect([1 1 1]);
            view(90,0)
            %axis tight
            axis off
            
            % filenames
            fname_fig = [txt_sym,'_Enorm_',kpt_txt,'_band_',num2str(nb),'.png'];
            pathFig = [datLoc,fBase,'\',fname_fig];
            saveas(gcf,pathFig);
                        
            % plot cross-sectional Ey profile
            EyExpr = 'emw.Ey';
            absEyExpr = ['abs(' EyExpr ')'];
            EyMax = mphmax(model,absEyExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',kn);
            ds.EyMax(kn,nb) = EyMax;
            %strplot.set('solnum', nb);
            Eyplt.set('looplevel',{num2str(nb) num2str(kn)});
            Eyplt_slc.set('expr', EyExpr);
            Eyplt_slc.set('rangedatamin', -EyMax).set('rangedatamax', EyMax);
            Eyplt_slc.set('rangecolormin', -EyMax).set('rangecolormax', EyMax);
            Eyplt_slc.set('colortable', 'WaveLight');
            Eyplt.run;
            
            figure;
            try
                mphplot(model,'Eyplt');
            catch err
            end
            title({['E-field, y component, ',txt_sol,', ',kpt_txt,'-point, '];...
                   ['band no. ',num2str(nb),', ',...
                   '\omega_o = ' num2str(ds.F(kn,nb)/1e12) ' THz']});
            colorbar('eastoutside');
            caxis([-EyMax EyMax]);
            colormap(colortable('WaveLight'));
            daspect([1 1 1]);
            view(90,0)
            %axis tight
            axis off
            
            % filenames
            fname_fig = [txt_sym,'_Ey_',kpt_txt,'_band_',num2str(nb),'.png'];
            pathFig = [datLoc,fBase,'\',fname_fig];
            saveas(gcf,pathFig);
            
            % plot cross-sectional Ez profile
            EzExpr = 'emw.Ez';
            absEzExpr = ['abs(' EzExpr ')'];
            EzMax = mphmax(model,absEzExpr,'volume','Dataset','pdset',...
                             'solnum',nb,'outersolnum',kn);
            ds.EzMax(kn,nb) = EzMax;
            %strplot.set('solnum', nb);
            Ezplt.set('looplevel',{num2str(nb) num2str(kn)});
            Ezplt_slc.set('expr', EzExpr);
            Ezplt_slc.set('rangedatamin', -EzMax).set('rangedatamax', EzMax);
            Ezplt_slc.set('rangecolormin', -EzMax).set('rangecolormax', EzMax);
            Ezplt_slc.set('colortable', 'WaveLight');
            Ezplt.run;
            
            figure;
            try
                mphplot(model,'Ezplt');
            catch err
            end
            title({['E-field, z component, ',txt_sol,', ',kpt_txt,'-point, '];...
                   ['band no. ',num2str(nb),', ',...
                   '\omega_o = ' num2str(ds.F(kn,nb)/1e12) ' THz']});
            colorbar('eastoutside');
            caxis([-EzMax EzMax]);
            colormap(colortable('WaveLight'));
            daspect([1 1 1]);
            view(90,0)
            %axis tight
            axis off
            
            % filenames
            fname_fig = [txt_sym,'_Ez_',kpt_txt,'_band_',num2str(nb),'.png'];
            pathFig = [datLoc,fBase,'\',fname_fig];
            saveas(gcf,pathFig);
            
            
            close all;
        end %of: for nb = 1:nbands
        %%
        % save fem data struct
        
        
        if ki == kpts %in: for nb = 1:nbands
            ds.femX = fem;
        end


    end %of: if P.saveplots && (ki==0 || ki==kpts) 
    
    

end %of: for ki = 0:kpts
%display('Postprocessing done');

if P.savedat
    mphfname = ['Bands_',txt_sym,'.mph'];
    path_mph = [datLoc,fBase,'\',mphfname];
    mphsave(model, path_mph);
end
ModelUtil.clear();
%display('Model and data saved. Run complete');
% rmpath(utilFilesPath)
end
