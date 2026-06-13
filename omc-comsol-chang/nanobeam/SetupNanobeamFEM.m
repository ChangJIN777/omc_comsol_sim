% Function to setup FEM simulation(s) to solve for optical and mechanical modes
%
% Inputs
% model: COMSOL model with geometry set up
% P: data structure with geometry and simulation params
% 
% Outputs
% model: updated COMSOL model with studies set up
% ds: data structure with FEM substructures

function [model,ds] = SetupNanobeamFEM(model,P)

% extract parameters from P
evenz = 0;
wid = P.w;
thi = P.th;
widMax = wid;

if P.solveMech
    % mechanical
    if isfield(P,'asymCav') && P.asymCav
        mevenx = 0;
        P.mevenx = mevenx;
    else
        mevenx = P.mevenx;
    end
    meveny = P.meveny;      % 1/-1 for symmetry/anti-symmetry; 0 for no symmetry
    mevenz = P.mevenz*strcmp(P.xsect,'rect');      % 1/-1 for symmetry/anti-symmetry; 0 for no symmetry
    P.mevenz = mevenz;
    
end
    
if P.solveOpt
    % optical
    if isfield(P,'asymCav') && P.asymCav
        oevenx = 0;
        P.oevenx = oevenx;
    end
    oevenx = P.oevenx^(P.holeatctr);      % 1/-1 for symmetry/anti-symmetry; 0 for no symmetry
    oeveny = P.oeveny;      % 1/-1 for symmetry/anti-symmetry; 0 for no symmetry
    oevenz = P.oevenz;      % 1/-1 for symmetry/anti-symmetry; 0 for no symmetry
    evenz = oevenz;
    airrad = P.airrad;      % radius of air cylinder around nanobeam
    nbeam = P.nbeam;        % refractive index of nanobeam
    lambda = P.lambda;      % target optical wavelength
    
    if airrad > widMax
        widMax = airrad;
    end
end

% extract geometry parameters for specifying boundary conditions

if isfield(P,'asymCav') && P.asymCav
    len = P.beamLen;    % for end of right half of beam
else
    len = P.beamLenHalf;        % length of beam
end
% totLen = len;

if P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML
    if isfield(P,'asymCav') && P.asymCav
        xL = -P.PMLLen;
        totLen = len + 2*P.PMLLen;
    else
        xL = 0;
        totLen = len + P.PMLLen;
    end
    
    if P.PMLLen > widMax
        widMax = P.PMLLen;
    end
end

% use half of beam height if z-symmetry enabled for rectangular nanobeam
% if strcmp(P.xsect,'rect') && abs(evenz)
%     thi = P.th/2;
% end

% geometry name
geomnames = fieldnames(mphmodel(model.geom));
geomname = geomnames{1};
beamgeom = model.geom(geomname);

%% Define material and properties
% define beam material
matTags = mphmodel(model.material);
if ~isfield(matTags,P.beamMat)
    bMat = model.material.create(P.beamMat);
else
    bMat = model.material(P.beamMat);
end
bMat.label(P.beamMat);
bMat_def = bMat.propertyGroup('def');

% mechanical properties
if P.solveMech

% density
mfem.rho = P.rho;
bMat_def.set('density', mfem.rho);

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
        error('ERROR: Insufficient number of stiffness constants specified');
    end
    bMat_aniso = bMat.propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
    bMat_aniso.set('DVo', P.D);
    
elseif (isfield(P,'E') && isfield(P,'nu')) 
    mfem.nu = P.nu;
    mfem.E = P.E;
    bMat_def.set('youngsmodulus', mfem.E);
    bMat_def.set('poissonsratio', mfem.nu);
else
    error('ERROR: No stiffness type specified');
end

% find index of domain corresponding to beam, then assign material to beam
% mfem.dia_domind = solid_index(beam,-1,[0,wid/2,thi/2],[0,0,1]); % beam
mfem.dia_domind = P.domSel.beam;
% mfem.pml_domind = [];
if isfield(P,'solveMechPML') && P.solveMechPML
    mfem.dia_domind = sort([mfem.dia_domind,P.domSel.PML]);
end
bMat.selection.geom(geomname, 3);    %to select domain
bMat.selection.set(mfem.dia_domind);
end % of if P.solveMech

% optical properties - conductivity, rel permittivity, refractive index
% first/second element of cell array refers to value for air/beam (e.g. diamond)
if P.solveOpt
% create air material node, if not already created
if ~isfield(matTags,'air')
    air = model.material.create('air');
else
    air = model.material('air');
end
air.label('Air');
air_def = air.propertyGroup('def');

ofem.sigma = {0,'1e-12[S/m]'};  % conductivity
ofem.epsilonr = {1,nbeam^2};    % relative permittivity
ofem.n = {1,nbeam};             % refractive index

% define optical properties of air and beam materials
air_ref = air.propertyGroup.create('air_ref', 'Refractive index');
air_def.set('electricconductivity', ofem.sigma{1});
air_def.set('relpermittivity', ofem.epsilonr{1});
air_def.set('relpermeability', 1);
air_ref.set('n', ofem.n{1});

bMat_ref = bMat.propertyGroup.create('bMat_ref', 'Refractive index');
bMat_def.set('electricconductivity', ofem.sigma{2});
bMat_def.set('relpermittivity', ofem.epsilonr{2});
bMat_def.set('relpermeability', 1);
bMat_ref.set('n', ofem.n{2});

% apply material to each domain
ofem.dia_domind = P.domSel.beam;
if ~P.solveMech
    bMat.selection.geom(geomname, 3);    %to select domain
    bMat.selection.set(ofem.dia_domind);
end

ofem.air_domind = P.domSel.cyl;
air.selection.geom(geomname, 3);    %to select domain
air.selection.set(ofem.air_domind);

% % pos parameter setting in solid_index:
% % air cylinder (nanobeam) is above (below) top boundary of nanobeam 
% % defined by p0 = [0 wid/2 thi] and n = [0 0 1]
% % ==> set pos = 1 (-1) for air cylinder (nanobeam).
% ofem.dia_domind = solid_index(beam,-1,[0,wid/2,thi/2],[0,0,1]);
% bMat.selection.geom(geomname, 3);    %to select domain
% bMat.selection.set(ofem.dia_domind);
% ofem.air_domind = solid_index(beam, 1,[0,wid/2,thi/2],[0,0,1]);
% air.selection.geom(geomname, 3);    %to select domain
% air.selection.set(ofem.air_domind);
end % of if P.solveOpt

display(['Beam material: ',P.beamMat,' added']);

%% Set up solid mechanics and define boundary conditions
if P.solveMech
% lemmx = linear elastic material, material no. x
% set material to Anisotropic, Voigt notation, and get elasticity and
% density from material definitions (above)
smech = model.physics.create('smech', 'SolidMechanics', geomname);
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
smech.selection.set(mfem.dia_domind);

% boundary conditions: all BCs set to free by default
% clear bnds


% fixed BCs - end of beam or PML pad
bnds.fixed_inds = [];
fixedBCs = smech.create('fixedBCs', 'Fixed', 2);
fixedBCs.label('Fixed Constraint');
if ~(isfield(P,'solveMechPML') && P.solveMechPML)
    bnds.fixed_inds = [bnds.fixed_inds,P.bndSel.beamXend];
end
if mevenx == 0
    if isfield(P,'solveMechPML') && P.solveMechPML
        % native PML: outer face is traction-free by default — no fixed BC needed
    else
        bnds.fixed_inds = [bnds.fixed_inds,P.bndSel.beamXsym];
    end
end
if (~isempty(bnds.fixed_inds))
    fixedBCs.selection.set(bnds.fixed_inds);
else
    fixedBCs.active(false);
end


% symmetric BCs
symBCs = smech.create('symBCs', 'SymmetrySolid', 2);
symBCs.label('Symmetric BC');
bnds.sym_inds = [];
if mevenx == 1
    bnds.sym_inds = [bnds.sym_inds,P.bndSel.beamXsym];
end
if meveny == 1
    bnds.sym_inds = [bnds.sym_inds,P.bndSel.beamYsym];
    if isfield(P,'solveMechPML') && P.solveMechPML
        bnds.sym_inds = [bnds.sym_inds,P.bndSel.PMLYsym];
    end
end
if mevenz == 1 && strcmp(P.xsect,'rect')
    bnds.sym_inds = [bnds.sym_inds,P.bndSel.beamZsym];
    if isfield(P,'solveMechPML') && P.solveMechPML
        bnds.sym_inds = [bnds.sym_inds,P.bndSel.PMLZsym];
    end
end
if (~isempty(bnds.sym_inds))
    symBCs.selection.set(bnds.sym_inds);
else
    symBCs.active(false);
end


% anti-symmetric BCs
asymBCs = smech.create('asymBCs', 'Antisymmetry', 2);
asymBCs.label('Anti-symmetric BC');
bnds.asym_inds = [];
if mevenx == -1
    bnds.asym_inds = [bnds.asym_inds,P.bndSel.beamXsym];
end
if meveny == -1
    bnds.asym_inds = [bnds.asym_inds,P.bndSel.beamYsym];
    if isfield(P,'solveMechPML') && P.solveMechPML
        bnds.asym_inds = [bnds.asym_inds,P.bndSel.PMLYsym];
    end
end
if mevenz == -1 && strcmp(P.xsect,'rect')
    bnds.asym_inds = [bnds.asym_inds,P.bndSel.beamZsym];
    if isfield(P,'solveMechPML') && P.solveMechPML
        bnds.asym_inds = [bnds.asym_inds,P.bndSel.PMLZsym];
    end
end
if (~isempty(bnds.asym_inds))
    asymBCs.selection.set(bnds.asym_inds);
else
    asymBCs.active(false);
end

mfem.bnds = bnds;

if isfield(P,'solveMechPML') && P.solveMechPML
    pml = model.coordSystem.create('pml1', geomname, 'PML');
    pml.selection.set(P.domSel.PML);
    pml.set('ScalingType', 'userDefined');
    pml.set('directions', '2');
    pml.setIndex('dmax', '1[mm]', 0);
    pml.setIndex('dmax', '1[mm]', 1);
    pml.set('wavelengthSourceType', 'userDefined');
    v_long = sqrt(P.D(1) / P.rho);
    lambda_mech = v_long / P.freq;
    pml.set('typicalWavelength', [num2str(lambda_mech), '[m]']);
end

display('Solid Mechanics added - boundary conditions done');
end % of if P.solveMech


%% Set up physics (RF EMW freq domain) and define boundary conditions
if P.solveOpt
emw = model.physics.create('emw', 'ElectromagneticWaves', geomname);
emw.selection.set([ofem.dia_domind, ofem.air_domind]);

% boundary conditions: default = PEC for boundaries adjacent to single domain
% - not applicable for boundaries adjacent to two domains (i.e. between nanobeam and air cylinder)
% (this is the continuity condition in COMSOL 3.5)
clear bnds

% Perfect magnetic conductor (PMC) - set for even modes
PMCbnd = emw.create('PMCbnd', 'PerfectMagneticConductor', 2);
PMCbnd.label('Perfect Magnetic Conductor');
bnds.pmc_inds = [];
if oevenx == 1
    %set BC for faces in yz-plane containing point (0,0,0)
    bnds.pmc_inds = [bnds.pmc_inds,P.bndSel.cylXsym];
end
if oeveny == 1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.pmc_inds = [bnds.pmc_inds,P.bndSel.cylYsym];
end
if oevenz == 1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.pmc_inds = [bnds.pmc_inds,P.bndSel.cylZsym];
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
bnds.scat_inds = [P.bndSel.cylXend,P.bndSel.cylCurv];
if oevenx == 0
    %set BC for faces at left end of beam
    bnds.scat_inds = [bnds.scat_inds,P.bndSel.cylXsym];
end
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
if oevenx == -1
    %set BC for faces in yz-plane containing point (0,0,0)
    bnds.pec_inds = [bnds.pec_inds,P.bndSel.cylXsym];
end
if oeveny == -1
    %set BC for faces in xz-plane containing point (0,0,0)
    bnds.pec_inds = [bnds.pec_inds,P.bndSel.cylYsym];
end
if oevenz == -1 && strcmp(P.xsect,'rect')
    %set BC for faces in xy-plane containing point (0,0,0)
    % for rectangular cross-section only
    bnds.pec_inds = [bnds.pec_inds,P.bndSel.cylZsym];
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

ofem.bnds = bnds;

display('EM Waves Frequency Domain added - boundary conditions done');
end % of if P.solveOpt


%% data structure management
ds.P = P;
if P.solveMech
    ds.mfem = mfem;
end
if P.solveOpt
    ds.ofem = ofem;
end
end
